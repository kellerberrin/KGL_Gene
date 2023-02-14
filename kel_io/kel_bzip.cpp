//
// Created by kellerberrin on 13/12/20.
//

#include "kel_bzip.h"
#include "kel_exec_env.h"
#include "kel_utility.h"

#include <zlib.h>
#include <fstream>


namespace kel = kellerberrin;

// Stops processing and empties queues, the BGZReader can again be re-opened after being closed.
void kel::BGZReader::close() {

  record_counter_ = 0;
  close_stream_ = true; // Stop processing and set eof condition.
  while (not readLine().EOFRecord()); // Drain the queues.
  bgz_file_.close(); // Close the physical file.
  reader_state_ = BGZReaderState::STOPPED;

}

// Begin processing and decompressing.
bool kel::BGZReader::open(const std::string &file_name) {

  // Cannot re-open an active object.
  if (reader_state_ == BGZReaderState::ACTIVE) {

    ExecEnv::log().error("BGZReader::open; stream is already active; call close().");
    return false;

  }

  record_counter_ = 0;
  file_name_ = file_name;
  close_stream_ = false;
  line_eof_ = false;

  try {

    // Open input '.bgz' file in a try block to catch any unexpected I/O exceptions.
    bgz_file_.open(file_name_, std::ios::binary | std::ios::ate);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::open; I/O error; could not open file: {}", file_name);
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("BGZReader::open; File: {} unexpected I/O exception: {}", file_name, e.what());
    return false;

  }

  // File is open so start processing.
  // Start decompressing bgz blocks.
  reader_thread_.enqueueVoid(&BGZReader::decompressGZBlockFile, this);
  // Unpacks the decompressed data and queues line records.
  line_assemble_thread_.enqueueVoid(&BGZReader::assembleRecords, this);

  return true;

}


std::optional<std::unique_ptr<kel::BaseStreamIO>> kel::BGZReader::getStreamIO( const std::string& file_name
                                                                             , size_t decompression_threads) {

  auto stream_ptr = std::make_unique<BGZReader>(decompression_threads);
  if (stream_ptr->open(file_name)) {

    return stream_ptr;

  }

  ExecEnv::log().error("BGZReader::getStreamIO; error opening file: {}", file_name);
  return std::nullopt;

}


// Remove decompressed line records from the BGZreader.
kel::IOLineRecord kel::BGZReader::readLine() {

// Dont block if eof reached.
  if (line_eof_ and line_queue_.empty()) return IOLineRecord::createEOFMarker();
// Return next available sequential line record.
  return line_queue_.waitAndPop();

}

// Physically read the '.bgz' file, check the file structure for errors, and submit compressed data blocks to be decompressed.
void kel::BGZReader::decompressGZBlockFile() {

  if (not bgz_file_.good()) {

    ExecEnv::log().error("BGZReader::readDecompressFile; failed to open bgz file: {}", file_name_);
    return;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file_.tellg();
  bgz_file_.seekg(0);

  size_t file_offset{0};
  size_t block_count{0};

  while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_) and bgz_file_.good()) {

    // Close down the stream gracefully.
    if (close_stream_) {

      // Push the eof marker.
      decompress_queue_.push(decompress_threads_.enqueueFuture(&BGZReader::decompressBlock, nullptr));
      return;

    }

    ++block_count;

    // Nasty but necessary.
    //  Read the header block.
    auto read_vector_ptr = std::make_shared<CompressedBlock>();
    bgz_file_.read(reinterpret_cast<char *>(&(read_vector_ptr->header_block_)), HEADER_SIZE_);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::readDecompressFile; Block {}, header file read error", block_count);
      // nullptr signals EOF downstream.
      decompress_queue_.push(decompress_threads_.enqueueFuture(&BGZReader::decompressBlock, nullptr));
      return;

    }

    // Nasty but necessary.
    // Read the compressed data
    size_t compressed_data_size = read_vector_ptr->header_block_.block_size + 1;
    char *read_ptr = reinterpret_cast<char *>(&(read_vector_ptr->compressed_block_[HEADER_SIZE_]));
    size_t read_data_size = compressed_data_size - HEADER_SIZE_;

    // Nasty but necessary.
    bgz_file_.read(read_ptr, read_data_size);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::readDecompressFile; Block {}, data file read error", block_count);
      // nullptr signals EOF downstream.
      decompress_queue_.push(decompress_threads_.enqueueFuture(&BGZReader::decompressBlock, nullptr));
      return;

    }

    read_vector_ptr->block_id_ = block_count;
    read_vector_ptr->data_size_ = compressed_data_size;

    file_offset += compressed_data_size;

    // Decompress the data.
    std::future<DecompressedType> future = decompress_threads_.enqueueFuture(&BGZReader::decompressBlock, read_vector_ptr);
    // Queue for processing.
    decompress_queue_.push(std::move(future));

  } // while blocks are available

  size_t remaining_chars = bgz_file_size - bgz_file_.tellg();
  checkEOFMarker(remaining_chars);

  // nullptr signals EOF downstream.
  decompress_queue_.push(decompress_threads_.enqueueFuture(&BGZReader::decompressBlock, nullptr));

}

void kel::BGZReader::checkEOFMarker(size_t remaining_chars) {

  // Not terminated so verify the EOF block.
  if (remaining_chars != EOF_MARKER_SIZE_) {

    ExecEnv::log().error("BGZReader::readDecompressFile; EOF Remaining bytes: {}, expected EOF remaining bytes: {}"
                         , remaining_chars, EOF_MARKER_SIZE_);
  } else {

    uint8_t eof_marker[EOF_MARKER_SIZE_];
    bgz_file_.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

    for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

      if (EOF_MARKER_[index] != eof_marker[index]) {

        ExecEnv::log().error("BGZReader::readDecompressFile; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                             index, eof_marker[index], EOF_MARKER_[index]);
      }

    }

  }

}


kel::BGZReader::DecompressedType  kel::BGZReader::decompressBlock(std::shared_ptr<CompressedBlock> compressed_data) {

  // Check for EOF and pass onto assembleRecords().
  if (not compressed_data) {

    return nullptr;

  }

  // Not EOF so decompress the block.
  z_stream_s zlib_params;
  char uncompressed_data[MAX_UNCOMPRESSED_SIZE_];
  size_t uncompressed_size{0};

  // Inflate the compressed data.
  zlib_params.next_in = reinterpret_cast<unsigned char*>(&(compressed_data->compressed_block_));
  zlib_params.avail_in = compressed_data->data_size_;
  zlib_params.next_out = reinterpret_cast<unsigned char*>(&uncompressed_data[0]);
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
    ExecEnv::log().error("BGZReader::decompressBlock; ::inflateInit2() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    return nullptr;
  }

  return_code = ::inflate(&zlib_params, Z_FINISH);
  if (return_code != Z_STREAM_END) {

    std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
    ExecEnv::log().error("BGZReader::decompressBlock; ::inflate() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    return nullptr;
  }

  return_code = ::inflateEnd(&zlib_params);
  if (return_code < Z_OK) {

    std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
    ExecEnv::log().error("BGZReader::decompressBlock; ::inflateEnd() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                         return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
    return nullptr;
  }

  // Parse the uncompressed block as much as possible.
  // The leading and trailing partial records are processed in the next stage (if not complete).

  uncompressed_size = MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out;

  auto uncompressed_ptr = std::make_unique<DecompressedBlock>();
  uncompressed_ptr->block_id = compressed_data->block_id_;

  std::string_view block_view(uncompressed_data, uncompressed_size);
  // Important note; if a final '\n' then the view parser allocates an empty string view.
  std::vector<std::string_view> view_vector = Utility::viewTokenizer(block_view, EOL_MARKER_);
  for (auto const& view : view_vector) {

    uncompressed_ptr->parsed_records.push_back(std::make_unique<std::string>(view));

  }

  return uncompressed_ptr;

}

void kel::BGZReader::assembleRecords() {

  record_counter_ = 0;
  size_t block_count{0};
  std::unique_ptr<std::string> previous_line_record;

  while(true) {

    std::future<DecompressedType> future = decompress_queue_.waitAndPop();
    DecompressedType block = future.get();
    // Check for EOF.
    if (not block) {
      // Push the EOF marker.
      line_queue_.push(IOLineRecord::createEOFMarker());
      line_eof_ = true;
      return;
    }

    // Check the decompressed blocks are sequential.
    ++block_count;
    if (block->block_id != block_count) {

      ExecEnv::log().warn("BGZReader::assembleRecords; Block mismatch, Queued block: {}, Counted block: {}", block->block_id, block_count);

    }

      // If the previous line record is defined, then concatenate with the first line record and queue.
    if (previous_line_record and not block->parsed_records.empty()) {

      previous_line_record->append(*(block->parsed_records.front()));

      if ( block->parsed_records.size() >= 2) {

        ++record_counter_;
        line_queue_.push(IOLineRecord(record_counter_, std::move(*previous_line_record)));
        previous_line_record = nullptr;

      }

    } else if (not block->parsed_records.empty()) {

      ++record_counter_;
      line_queue_.push(IOLineRecord(record_counter_, std::move(*(block->parsed_records.front()))));
      previous_line_record = nullptr;

    }

    // Queue the 2nd to 2nd last records.
    size_t line_count = block->parsed_records.size();
    for (size_t index = 1; index < line_count-1; ++index) {

      ++record_counter_;
      line_queue_.push(IOLineRecord(record_counter_, std::move(*(block->parsed_records[index]))));

    }

    if (not previous_line_record and not block->parsed_records.empty()) {

      // If this block is complete then queue the line and clear the previous_line_record.
      // Store the last record in the previous_line_record.
      previous_line_record = std::move(block->parsed_records.back());

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



