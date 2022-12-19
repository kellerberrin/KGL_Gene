//
// Created by kellerberrin on 13/12/20.
//

#include "zlib.h"

#include <fstream>

#include "kel_bzip.h"
#include "kel_exec_env.h"
#include "kel_utility.h"

namespace kel = kellerberrin;


bool kel::BGZReader::close() {

  shutdown_ = true;
  // Empty the queues to flush through the eof markers and shutdown the worker threads.
  if (not line_queue_.empty()) while(readLine());
  // The std::future may not be active, so eat any exceptions thrown by querying it.
  try {

    reader_return_.get();

  } catch(...) {}

  bgz_file_.close();

  // Reset the flags
  shutdown_ = false;
  line_eof_ = false;
  decompression_error_ = false;

  return true;

}


bool kel::BGZReader::open(const std::string &file_name) {

  try {

    record_counter_ = 0;
    file_name_ = file_name;
    // Open input file.

    bgz_file_.open(file_name_, std::ios::binary | std::ios::ate);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::open; I/O error; could not open file: {}", file_name);
      decompression_error_ = true;
      return false;

    }
  }
  catch (std::exception const &e) {

    ExecEnv::log().error("BGZReader::open; File: {} unexpected I/O exception: {}", file_name, e.what());
    decompression_error_ = true;
    return false;

  }

  // Reset the flags
  shutdown_ = false;
  line_eof_ = false;
  decompression_error_ = false;
  // File is open so start processing.
  reader_return_ = reader_thread_.enqueueTask(&BGZReader::decompressGZBlockFile, this);

  return true;

}


kel::IOLineRecord kel::BGZReader::readLine() {

// Dont block if eof reached.
  if (line_eof_ and line_queue_.empty()) return QUEUED_EOF_MARKER;
// Return next available sequential line record.
  return line_queue_.waitAndPop();

}


bool kel::BGZReader::verify(const std::string &file_name, bool silent) {

  std::ifstream bgz_file(file_name, std::ios::binary | std::ios::ate);

  if (not bgz_file.good()) {

    ExecEnv::log().warn("BGZReader::verify; Problem opening file: {}", file_name);
    return false;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file.tellg();
  ExecEnv::log().info("Verifing bgz file: {}:, Size: {}", file_name, bgz_file_size);
  // Reset the file stream pointer to the file begining.
  bgz_file.seekg(0);

  size_t file_offset{0};
  size_t block_count{0};
  size_t total_compressed_size{0};
  size_t total_uncompressed_size{0};
  while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_) and not bgz_file.eof()) {

    ++block_count;
    // Read header.
    GZHeaderblock header_block;
    bgz_file.read(reinterpret_cast<char*>(&header_block), HEADER_SIZE_);
    file_offset += HEADER_SIZE_;
    // Check the header values.
    if (header_block.block_id_1 != BLOCK_ID1_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, bad header block id 1: {}, expected: {}",
                                            block_count, header_block.block_id_1, BLOCK_ID1_);
      return false;

    }

    if (header_block.block_id_2 != BLOCK_ID2_) {

      if (not silent) ExecEnv::log().error("BGZReader::verify; Block count: {}, bad header block id 2: {}, expected: {}",
                                           block_count, header_block.block_id_2, BLOCK_ID2_);
      return false;

    }

    if (header_block.subfield_id_1 != SUBFIELD_ID1_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                                            block_count, header_block.subfield_id_1, SUBFIELD_ID1_);
      return false;

    }

    if (header_block.subfield_id_2 != SUBFIELD_ID2_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                                            block_count, header_block.subfield_id_2, SUBFIELD_ID2_);
      return false;

    }

    if (header_block.length_extra_blocks != EXTRA_LENGTH_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, extra block size: {}, expected: {}",
                                            block_count, header_block.length_extra_blocks, EXTRA_LENGTH_);
      return false;

    }

    // skip the compressed data
    size_t compressed_data_size = header_block.block_size - BLOCK_SIZE_ADJUST_;
    if (compressed_data_size > MAX_UNCOMPRESSED_SIZE_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, Compressed block size: {} exceeds max Ccompressed size: {}",
                                            block_count, compressed_data_size, MAX_UNCOMPRESSED_SIZE_);
      return false;

    }
    bgz_file.seekg(compressed_data_size , std::ios_base::cur);
    file_offset += compressed_data_size;
    total_compressed_size += compressed_data_size;

    // Read the trailer block
    GZTrailerBlock trailer_block;
    bgz_file.read(reinterpret_cast<char*>(&trailer_block), TRAILER_SIZE_);
    file_offset += TRAILER_SIZE_;
    total_uncompressed_size += trailer_block.uncompressed_size;

    if (trailer_block.uncompressed_size > MAX_UNCOMPRESSED_SIZE_) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, Uncompressed block size: {} exceeds max uncompressed size: {}",
                                            block_count, trailer_block.uncompressed_size, MAX_UNCOMPRESSED_SIZE_);
      return false;

    }

    size_t file_position = static_cast<size_t>(bgz_file.tellg());
    if (file_position != file_offset) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; Block count: {}, file tellg: {}, calc file_offset: {}",
                                            block_count, bgz_file.tellg(), file_offset);
      return false;

    }

  } // while.


  // Check EOF string.
  size_t remaining_chars = bgz_file_size - bgz_file.tellg();
  if (remaining_chars != EOF_MARKER_SIZE_) {

    if (not silent) ExecEnv::log().error( "BGZReader::verify; Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                                          block_count, remaining_chars, EOF_MARKER_SIZE_);
    return false;
  }

  uint8_t eof_marker[EOF_MARKER_SIZE_];
  bgz_file.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

  for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

    if (EOF_MARKER_[index] != eof_marker[index]) {

      if (not silent) ExecEnv::log().error( "BGZReader::verify; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                                            index, eof_marker[index], EOF_MARKER_[index]);
      return false;

    }

  }

  ExecEnv::log().info("Verified: {}, Blocks: {}, Data Uncompressed: {}, Compressed: {}",
                      file_name,  block_count, total_uncompressed_size, total_compressed_size);

  return true;

}


bool kel::BGZReader::decompressGZBlockFile() {

  WorkflowThreads thread_pool(thread_count_);
  WorkflowThreads line_thread(1);

  // Unpacks the decompressed data and queues line records.
  line_thread.enqueueWork(&BGZReader::assembleRecords, this);

  if (not bgz_file_.good()) {

    ExecEnv::log().error("BGZReader::decompressGZBlockFile; failed to open bgz file: {}", file_name_);
    decompression_error_ = true;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file_.tellg();
  bgz_file_.seekg(0);

  size_t file_offset{0};
  size_t block_count{0};
  size_t total_compressed_size{0};

  while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_)
         and not bgz_file_.eof()
         and not shutdown_
         and not decompression_error_) {

    ++block_count;

    //  Must be created for each iteration
    std::shared_ptr<std::vector<std::byte>> read_vector_ptr(std::make_shared<std::vector<std::byte>>(MAX_UNCOMPRESSED_SIZE_));
    bgz_file_.read(reinterpret_cast<char *>(&(read_vector_ptr->at(0))), HEADER_SIZE_);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::decompressGZBlockFile; Block {}, header file read error", block_count);
      decompression_error_ = true;
      break;

    }

    // Nasty but necessary.
    auto header_block_ptr = reinterpret_cast<GZHeaderblock *>(&(read_vector_ptr->at(0)));
    // Read the compressed data
    size_t compressed_data_size = header_block_ptr->block_size + 1;
    char *read_ptr = reinterpret_cast<char *>(&(read_vector_ptr->at(HEADER_SIZE_)));
    size_t read_data_size = compressed_data_size - HEADER_SIZE_;

    // Nasty but necessary.
    bgz_file_.read(read_ptr, read_data_size);
    if (not bgz_file_.good()) {

      ExecEnv::log().error("BGZReader::decompressGZBlockFile; Block {}, data file read error", block_count);
      decompression_error_ = true;
      break;

    }

    file_offset += compressed_data_size;
    total_compressed_size += compressed_data_size;


    // Decompress the data.
    std::future<UncompressedBlock> future = thread_pool.enqueueTask(&BGZReader::decompressBlock,
                                                                    this,
                                                                    block_count,
                                                                    read_vector_ptr,
                                                                    compressed_data_size,
                                                                    false);
    if (decompression_error_) {

      break;

    }
    // Queue for processing.
    decompress_queue_.push(std::move(future));

  } // while blocks available and not shutdown.

  size_t remaining_chars = bgz_file_size - bgz_file_.tellg();

  // If terminated just return.
  if (shutdown_ or decompression_error_) {

    // Queue an empty block to signal an EOF or error condition.
    std::future<UncompressedBlock> future = thread_pool.enqueueTask(&BGZReader::decompressBlock, this, 0, nullptr, 0, true);
    // Queue for processing.
    decompress_queue_.push(std::move(future));

    return true;

  }

  // Not terminated so verify the EOF block.
  if (remaining_chars != EOF_MARKER_SIZE_) {

    ExecEnv::log().error("BGZReader::decompressGZBlockFile; Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                         block_count, remaining_chars, EOF_MARKER_SIZE_);
  } else {

    uint8_t eof_marker[EOF_MARKER_SIZE_];
    bgz_file_.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

    for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

      if (EOF_MARKER_[index] != eof_marker[index]) {

        ExecEnv::log().error("BGZReader::decompressGZBlockFile; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                             index, eof_marker[index], EOF_MARKER_[index]);
      }

    }

  }

  // Queue an empty block to signal an EOF or error condition.
  std::future<UncompressedBlock> future = thread_pool.enqueueTask(&BGZReader::decompressBlock, this, 0, nullptr, 0, true);
  // Queue for processing.
  decompress_queue_.push(std::move(future));

  return true;

}



kel::UncompressedBlock kel::BGZReader::decompressBlock(size_t block_count,
                                                       std::shared_ptr<std::vector<std::byte>> compressed_data,
                                                       size_t compressed_data_size,
                                                       bool eof_flag) {

  UncompressedBlock uncompressed_block;
  z_stream_s zlib_params;
  char uncompressed_data[MAX_UNCOMPRESSED_SIZE_];
  size_t uncompressed_size{0};

  // Decompress or otherwise return an empty block (eof marker).
  if (not eof_flag) {

    // Inflate the compressed data.
    zlib_params.next_in = reinterpret_cast<unsigned char*>(&(compressed_data->at(0)));
    zlib_params.avail_in = compressed_data_size;
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
      decompression_error_ = true;

    }

    return_code = ::inflate(&zlib_params, Z_FINISH);
    if (return_code != Z_STREAM_END) {

      std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
      ExecEnv::log().error("BGZReader::decompressBlock; ::inflate() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                           return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
      decompression_error_ = true;

    }

    return_code = ::inflateEnd(&zlib_params);
    if (return_code < Z_OK) {

      std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
      ExecEnv::log().error("BGZReader::decompressBlock; ::inflateEnd() fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                           return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);
      decompression_error_ = true;

    }

    // Parse the uncompressed block as much as possible.
    // The leading and trailing partial records are processed in the next stage (if not complete).
    if (not decompression_error_) {

      uncompressed_size = MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out;
      uncompressed_block.block_id = block_count;
      uncompressed_block.eof_flag = false;

      std::string_view block_view(uncompressed_data, uncompressed_size);
      // Important note; if a final '\n' then the view parser allocates an empty string view.
      std::vector<std::string_view> view_vector = Utility::view_tokenizer(block_view, EOL_MARKER_);
      for (auto const& view : view_vector) {

        uncompressed_block.parsed_records.push_back(std::make_unique<std::string>(view));

      }

    } else {

      uncompressed_block.eof_flag = true;

    }

  } else {

    uncompressed_block.eof_flag = true;

  }


  return uncompressed_block;

}


void kel::BGZReader::assembleRecords() {

  record_counter_ = 0;
  size_t block_count{0};
  std::unique_ptr<std::string> previous_line_record;

  while(true) {

    std::future<UncompressedBlock> future = decompress_queue_.waitAndPop();
    UncompressedBlock block = future.get();
    // Check for eof.
    if (block.eof_flag) {

      line_eof_ = true;
      break;

    }

    ++block_count;

    if (block.block_id != block_count) {

      ExecEnv::log().warn("BGZReader::assembleRecords; Block mismatch, Queued block: {}, Counted block: {}", block.block_id, block_count);

    }

      // If the previous line record is defined, then concatenate with the first line record and queue.
    if (previous_line_record and not block.parsed_records.empty()) {

      previous_line_record->append(*block.parsed_records.front());

      if ( block.parsed_records.size() >= 2) {

        ++record_counter_;
        line_queue_.push(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(previous_line_record)));
        previous_line_record = nullptr;

      }

    } else if (not block.parsed_records.empty()) {

      ++record_counter_;
      line_queue_.push(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(block.parsed_records.front())));
      previous_line_record = nullptr;

    }

    // Queue the 2nd to 2nd last last records.
    size_t line_count = block.parsed_records.size();
    for (size_t index = 1; index < line_count-1; ++index) {

      ++record_counter_;
      line_queue_.push(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(block.parsed_records[index])));

    }

    if (not previous_line_record and not block.parsed_records.empty()) {

      // If this block is complete then queue the line and clear the previous_line_record.
      // Store the last record in the previous_line_record.
      previous_line_record = std::move(block.parsed_records.back());

    }

  } // while.

  // Queue the final record if found and non-empty
  if (previous_line_record) {

    if (not previous_line_record->empty()) {

      ++record_counter_;
      line_queue_.push(std::pair<size_t, std::unique_ptr<std::string>>(record_counter_, std::move(previous_line_record)));

    }

  }

  // Push the eof marker.
  line_queue_.push(QUEUED_EOF_MARKER);

}



