//
// Created by kellerberrin on 17/01/23.
//

#include "kel_bzip_workflow.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <fstream>
#include <zlib.h>


namespace kel = kellerberrin;


bool kel::BGZStream::close() {

  close_stream_ = true; // stop processing
  while (not readLine().EOFRecord()); // drain the queues.
  bgz_file_.close();
  stream_state_ = BGZStreamState::STOPPED;
  return true;

}


bool kel::BGZStream::open(const std::string &file_name) {

  // Cannot re-open the object.
  if (stream_state_ == BGZStreamState::ACTIVE) {

    ExecEnv::log().error("BGZStream::open; stream is already active; call close().");
    return false;

  }

  record_counter_ = 0;
  file_name_ = file_name;
  close_stream_ = false;
  line_eof_ = false;

  try {

    // Open input file.
    bgz_file_.open(file_name_, std::ios::binary | std::ios::ate);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZStream::open; I/O error; could not open file: {}", file_name);
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("BGZStream::open; File: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  // File is open so start processing.

  // Activate the decompression workflow.
  decompression_workflow_.activateWorkflow(decompression_threads_, &BGZStream::decompressBlock, this);
  // Begin decompressing blocks of data.
  reader_return_ = reader_thread_.enqueueFuture(&BGZStream::readDecompressFile, this);
  // Begin queueing decompressed text records.
  assemble_records_thread_.enqueueVoid(&BGZStream::assembleRecords, this);
  // Set the object state.
  stream_state_ = BGZStreamState::ACTIVE;

  return true;

}


kel::IOLineRecord kel::BGZStream::readLine() {

  // Dont block if eof reached.
  if (line_eof_ and line_queue_.empty()) return IOLineRecord::createEOFMarker();
  // Return next available sequential line record.
  return line_queue_.waitAndPop();

}


bool kel::BGZStream::readDecompressFile() {

  if (not bgz_file_.good()) {

    ExecEnv::log().error("BGZStream::readDecompressFile; failed to open bgz file: {}", file_name_);
    return false;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file_.tellg();
  bgz_file_.seekg(0);

  size_t file_offset{0};
  size_t block_count{0};

  while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_) and not bgz_file_.eof()) {

    // Close down the stream gracefully.
    if (close_stream_) {

      // Push the eof marker.
      decompression_workflow_.push(nullptr);
      return true;

    }

    ++block_count;
    auto compressed_ptr = readCompressedBlock(block_count);
    if (not compressed_ptr->io_success_) {

      ExecEnv::log().error("BGZStream::readDecompressFile; Encountered I/O error reading compressed block {}", block_count);
      return false;
    }

    file_offset += compressed_ptr->data_size_;
    decompression_workflow_.push(std::move(compressed_ptr));

  } // While compressed blocks available

  // Push the workflow stop token.
  decompression_workflow_.push(nullptr);
  // Verify the trailing block.
  size_t remaining_chars = bgz_file_size - bgz_file_.tellg();
  return checkEOFMarker(remaining_chars);

}


bool kel::BGZStream::checkEOFMarker(size_t remaining_chars) {

  // Not terminated so verify the EOF block.
  if (remaining_chars != EOF_MARKER_SIZE_) {

    ExecEnv::log().error("BGZStream::readDecompressFile; EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                         remaining_chars, EOF_MARKER_SIZE_);
    return false;

  } else {

    uint8_t eof_marker[EOF_MARKER_SIZE_];
    bgz_file_.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

    for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

      if (EOF_MARKER_[index] != eof_marker[index]) {

        ExecEnv::log().error("BGZStream::readDecompressFile; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                             index, eof_marker[index], EOF_MARKER_[index]);
        return false;

      }

    }

  }

  return true;

}



kel::BGZStream::CompressedType kel::BGZStream::readCompressedBlock(size_t block_count) {

  // First read the header.
  //  Must be created for each iteration
  std::unique_ptr<CompressedBlock> read_vector_ptr(std::make_unique<CompressedBlock>());
  read_vector_ptr->block_id_ = block_count;
  // Nasty but necessary.
  bgz_file_.read(reinterpret_cast<char *>(&(read_vector_ptr->header_block_)), HEADER_SIZE_);
  if (not bgz_file_.good()) {

    ExecEnv::log().error("BGZStream::readCompressedBlock; Block {}, header file read error", block_count);
    read_vector_ptr->io_success_ = false;
    return read_vector_ptr;

  }

  // Using the header block_size, now read the compressed byte data.
  // Nasty but necessary.
  int32_t total_block_size = read_vector_ptr->header_block_.block_size + 1;
  int32_t compressed_data_size = total_block_size - HEADER_SIZE_;
  // Nasty but necessary.
  char *read_ptr = reinterpret_cast<char *>(&(read_vector_ptr->compressed_block_[HEADER_SIZE_]));
  bgz_file_.read(read_ptr, compressed_data_size);
  if (not bgz_file_.good()) {

    ExecEnv::log().error("BGZStream::readCompressedBlock; Block {}, data file read error", block_count);
    read_vector_ptr->io_success_ = false;
    return read_vector_ptr;

  }

  read_vector_ptr->data_size_ = total_block_size;
  read_vector_ptr->io_success_ = true;

  return read_vector_ptr;

}


std::optional<kel::BGZStream::DecompressedType> kel::BGZStream::decompressBlock(CompressedType compressed_ptr) {

  // Check if a stop token.
  if (not compressed_ptr) {

    // just return the nullptr;
    return nullptr;

  }

  auto decompressed_ptr = std::make_unique<DecompressedBlock>();
  z_stream_s zlib_params;

  // Inflate the compressed data.
  zlib_params.next_in = reinterpret_cast<unsigned char*>(&(compressed_ptr->compressed_block_[0]));
  zlib_params.avail_in = compressed_ptr->data_size_;
  zlib_params.next_out = reinterpret_cast<unsigned char*>(&(decompressed_ptr->decompressed_data_[0]));
  zlib_params.avail_out = MAX_UNCOMPRESSED_SIZE_;
  zlib_params.data_type = Z_TEXT;
  zlib_params.total_out = 0;
  zlib_params.total_in = 0;
  zlib_params.zalloc = Z_NULL;
  zlib_params.zfree = Z_NULL;
  zlib_params.opaque = nullptr;

  int return_code = ::inflateInit2(&zlib_params, INFLATE_WINDOW_FLAG_);
  if (return_code < Z_OK) {

    std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
    ExecEnv::log().error("BGZStream::decompressBlock; ::inflateInit2() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    decompressed_ptr->decompress_success_ = false;
    return decompressed_ptr;

  }

  return_code = ::inflate(&zlib_params, Z_FINISH);
  if (return_code != Z_STREAM_END) {

    std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
    ExecEnv::log().error("BGZStream::decompressBlock; ::inflate() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    decompressed_ptr->decompress_success_ = false;
    return decompressed_ptr;

  }

  return_code = ::inflateEnd(&zlib_params);
  if (return_code < Z_OK) {

    std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
    ExecEnv::log().error("BGZStream::decompressBlock; ::inflateEnd() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    decompressed_ptr->decompress_success_ = false;
    return decompressed_ptr;

  }

  decompressed_ptr->data_size_ = MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out;
  decompressed_ptr->block_id_ = compressed_ptr->block_id_;
  decompressed_ptr->decompress_success_ = true;

  std::string_view block_view(&(decompressed_ptr->decompressed_data_[0]), decompressed_ptr->data_size_);
  // Important note; if a final '\n' then the view parser allocates an empty string view.
  std::vector<std::string_view> view_vector = Utility::viewTokenizer(block_view, EOL_MARKER_);
  for (auto const& view : view_vector) {

    decompressed_ptr->parsed_lines_.push_back(std::make_unique<std::string>(view));

  }

  return decompressed_ptr;

}


void kel::BGZStream::assembleRecords() {

  record_counter_ = 0;
  size_t block_count{0};
  std::unique_ptr<std::string> previous_line_record;

  while(true) {

    DecompressedType block_ptr = decompression_workflow_.waitAndPop();
    // Check for eof.
    if (not block_ptr) {

      line_eof_ = true;
      break;

    }

    if (not block_ptr->decompress_success_) {

      ExecEnv::log().warn("BGZStream::assembleRecords; decompress error with block: {}", block_ptr->block_id_);
      return;

    }

    ++block_count;
    if (block_ptr->block_id_ != block_count) {

      ExecEnv::log().warn("BGZStream::assembleRecords; Block mismatch, Queued block: {}, Counted block: {}", block_ptr->block_id_, block_count);

    }

    // If the previous line record is defined, then concatenate with the first line record and queue.
    if (previous_line_record and not block_ptr->parsed_lines_.empty()) {

      previous_line_record->append(*(block_ptr->parsed_lines_.front()));

      if ( block_ptr->parsed_lines_.size() >= 2) {

        ++record_counter_;
        line_queue_.push(IOLineRecord(record_counter_, std::move(*previous_line_record)));
        previous_line_record = nullptr;

      }

    } else if (not block_ptr->parsed_lines_.empty()) {

      ++record_counter_;
      line_queue_.push(IOLineRecord(record_counter_, std::move(*(block_ptr->parsed_lines_.front()))));
      previous_line_record = nullptr;

    }

    // Queue the 2nd to 2nd last last records.
    size_t line_count = block_ptr->parsed_lines_.size();
    for (size_t index = 1; index < line_count-1; ++index) {

      ++record_counter_;
      line_queue_.push(IOLineRecord(record_counter_, std::move(*(block_ptr->parsed_lines_[index]))));

    }

    if (not previous_line_record and not block_ptr->parsed_lines_.empty()) {

      // If this block is complete then queue the line and clear the previous_line_record.
      // Store the last record in the previous_line_record.
      previous_line_record = std::move(block_ptr->parsed_lines_.back());

    }

  } // while.

  // Queue the final record if found and non-empty
  if (previous_line_record) {

    if (not previous_line_record->empty()) {

      ++record_counter_;
      line_queue_.push(IOLineRecord(record_counter_, std::move(*previous_line_record)));

    }

  }

  // Push the eof marker.
  line_queue_.push(IOLineRecord::createEOFMarker());

}

