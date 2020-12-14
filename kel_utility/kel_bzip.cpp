//
// Created by kellerberrin on 13/12/20.
//

#include "zlib.h"

#include <fstream>

#include "kel_bzip.h"
#include "kel_exec_env.h"
#include "kel_thread_pool.h"

namespace kel = kellerberrin;


bool kel::GZBlockDecompression::verifyGZBlockFile() {

  // Open for binary reading.
  std::ifstream bgz_file(file_name_, std::ios::binary | std::ios::ate);

  if (not bgz_file.good()) {

    ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; failed to open bgz file: {}", file_name_);
    return false;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file.tellg();
  ExecEnv::log().info("Verifing bgz file: {}:, Size: {}", file_name_, bgz_file_size);
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

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, bad header block id 1: {}, expected: {}",
                           block_count, header_block.block_id_1, BLOCK_ID1_);
      return false;

    }

    if (header_block.block_id_2 != BLOCK_ID2_) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, bad header block id 2: {}, expected: {}",
                           block_count, header_block.block_id_2, BLOCK_ID2_);
      return false;

    }

    if (header_block.subfield_id_1 != SUBFIELD_ID1_) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                           block_count, header_block.subfield_id_1, SUBFIELD_ID1_);
      return false;

    }

    if (header_block.subfield_id_2 != SUBFIELD_ID2_) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                           block_count, header_block.subfield_id_2, SUBFIELD_ID2_);
      return false;

    }

    if (header_block.length_extra_blocks != EXTRA_LENGTH_) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, extra block size: {}, expected: {}",
                           block_count, header_block.length_extra_blocks, EXTRA_LENGTH_);
      return false;

    }

    // skip the compressed data
    size_t compressed_data_size = header_block.block_size - BLOCK_SIZE_ADJUST_;
    bgz_file.seekg(compressed_data_size , std::ios_base::cur);
    file_offset += compressed_data_size;
    total_compressed_size += compressed_data_size;

    // Read the trailer block
    GZTrailerBlock trailer_block;
    bgz_file.read(reinterpret_cast<char*>(&trailer_block), TRAILER_SIZE_);
    file_offset += TRAILER_SIZE_;
    total_uncompressed_size += trailer_block.uncompressed_size;

    if (trailer_block.uncompressed_size > MAX_UNCOMPRESSED_SIZE_) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, Uncompressed block size: {} exceeds max uncompressed size: {}",
                           block_count, trailer_block.uncompressed_size, MAX_UNCOMPRESSED_SIZE_);
      return false;

    }

    size_t file_position = static_cast<size_t>(bgz_file.tellg());
    if (file_position != file_offset) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Block count: {}, file tellg: {}, calc file_offset: {}",
                           block_count, bgz_file.tellg(), file_offset);
      return false;

    }

  } // while.

  // Check EOF string.
  size_t remaining_chars = bgz_file_size - bgz_file.tellg();
  if (remaining_chars != EOF_MARKER_SIZE_) {

    ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                        block_count, remaining_chars, EOF_MARKER_SIZE_);
    return false;
  }

  uint8_t eof_marker[EOF_MARKER_SIZE_];
  bgz_file.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

  for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

    if (EOF_MARKER_[index] != eof_marker[index]) {

      ExecEnv::log().error("GZBlockDecompression::verifyGZBlockFile; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                           index, eof_marker[index], EOF_MARKER_[index]);
      return false;

    }

  }

  ExecEnv::log().info("Successfully verified: {}, Blocks: {}, Uncompressed size: {}, Compressed size: {}",
                      file_name_,  block_count, total_uncompressed_size, total_compressed_size);

  return true;

}





bool kel::GZBlockDecompression::decompressGZBlockFile(size_t thread_count) {

  ThreadPool thread_pool(thread_count);
  ThreadPool line_thread(1);

  // Unpacks the decompressed data and queues line records.
  line_thread.enqueueWork(&GZBlockDecompression::queueLines, this);
  // Open for binary reading. Reposition the stream pointer to the end of file.
  std::ifstream bgz_file(file_name_, std::ios::binary | std::ios::ate);

  if (not bgz_file.good()) {

    ExecEnv::log().error("GZBlockDecompression::decompressGZBlockFile; failed to open bgz file: {}", file_name_);
    return false;

  }

  // Read file size.
  const size_t bgz_file_size = bgz_file.tellg();
  bgz_file.seekg(0);

  size_t file_offset{0};
  size_t block_count{0};
  size_t total_compressed_size{0};
  size_t total_uncompressed_size{0};

  while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_) and not bgz_file.eof()) {

    ++block_count;

    // Read header.
    auto read_buffer = new std::byte[MAX_UNCOMPRESSED_SIZE_];
    bgz_file.read(reinterpret_cast<char*>(&read_buffer[0]), HEADER_SIZE_);

    // Nasty but necessary.
    auto header_block_ptr = reinterpret_cast<GZHeaderblock*>(&read_buffer[0]);
    // Read the compressed data
    size_t compressed_data_size = header_block_ptr->block_size + 1;
    char* read_ptr = reinterpret_cast<char*>(&read_buffer[HEADER_SIZE_]);
    size_t read_data_size = compressed_data_size - HEADER_SIZE_;

    // Nasty but necessary.
    bgz_file.read(read_ptr, read_data_size);
    file_offset += compressed_data_size;
    total_compressed_size += compressed_data_size;

    // Decompress the data.
    std::future<UncompressedBlock> future = thread_pool.enqueueTask( &GZBlockDecompression::decompressBlock,
                                                                     block_count,
                                                                     read_buffer,
                                                                     compressed_data_size);
    // Queue for removal.
    decompress_queue_.push(std::move(future));


  } // while blocks available.

  // Check EOF string.
  size_t remaining_chars = bgz_file_size - bgz_file.tellg();
  if (remaining_chars != EOF_MARKER_SIZE_) {

    ExecEnv::log().error("GZBlockDecompression::decompressGZBlockFile; Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                         block_count, remaining_chars, EOF_MARKER_SIZE_);
  }

  uint8_t eof_marker[EOF_MARKER_SIZE_];
  bgz_file.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

  for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

    if (EOF_MARKER_[index] != eof_marker[index]) {

      ExecEnv::log().error("GZBlockDecompression::decompressGZBlockFile; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                           index, eof_marker[index], EOF_MARKER_[index]);
    }

  }

  // Queue an empty block to signal an EOF or error condition.
  ExecEnv::log().info("Enqueueing EOF");
  std::future<UncompressedBlock> future = thread_pool.enqueueTask( &GZBlockDecompression::decompressBlock,
                                                                   0,
                                                                   nullptr,
                                                                   0);
  // Queue for removal.
  decompress_queue_.push(std::move(future));

  ExecEnv::log().info("Processed: {}, Blocks: {}, Uncompressed size: {}, Compressed size: {}",
                      file_name_,  block_count, total_uncompressed_size, total_compressed_size);


  return true;

}



kel::UncompressedBlock kel::GZBlockDecompression::decompressBlock( size_t block_count,
                                                                   std::byte* compressed_data,
                                                                   size_t compressed_data_size) {

  UncompressedBlock uncompressed_block;
  z_stream_s zlib_params;

  // Decompress otherwise return an empty block (eof marker).
  if (compressed_data_size != 0 and compressed_data != nullptr) {

    uncompressed_block.uncompressed_data = std::make_unique<char[]>(MAX_UNCOMPRESSED_SIZE_);
    uncompressed_block.uncompressed_size = MAX_UNCOMPRESSED_SIZE_;
    // Inflate the compressed data.
    zlib_params.next_in = reinterpret_cast<unsigned char*>(compressed_data);
    zlib_params.avail_in = compressed_data_size;
    zlib_params.next_out = reinterpret_cast<unsigned char*>(uncompressed_block.uncompressed_data.get());
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
      ExecEnv::log().error("GZBlockDecompression::decompressBlock; fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                           return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);

    }

    return_code = ::inflate(&zlib_params, Z_FINISH);
    if (return_code != Z_STREAM_END) {

      std::string zlib_msg = zlib_params.msg != nullptr ? zlib_params.msg : "no msg";
      ExecEnv::log().error("GZBlockDecompression::decompressBlock; fail, return code: {}, msg: {}, Uncompressed: {}, Consumed: {}",
                           return_code, zlib_msg, MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out, zlib_params.total_in);

    }

    uncompressed_block.uncompressed_size = MAX_UNCOMPRESSED_SIZE_ - zlib_params.avail_out;
    uncompressed_block.block_id = block_count;

  }

  delete[] compressed_data;

  return uncompressed_block;

}


void kel::GZBlockDecompression::queueLines() {

  bool eof{false};
  size_t block_count{0};
  do {

    ++block_count;
    std::future<UncompressedBlock> future = decompress_queue_.waitAndPop();
    UncompressedBlock block = future.get();
    if (block.uncompressed_size == 0 and not block.uncompressed_data) {

      eof = true;

    } else {

      if (block.block_id != block_count) {

        ExecEnv::log().warn("GZBlockDecompression::queueLines; Block mismatch, Queued block: {}, Counted block: {}", block.block_id, block_count);

      }

      const char* buffer = block.uncompressed_data.get();
      std::vector<size_t> line_offsets;
      for (size_t index = 0; index < block.uncompressed_size; ++index) {

        if (buffer[index] == '\n') {

          line_offsets.push_back(index);

        }

      }

      size_t remain_alpha{0};
      if (not line_offsets.empty()) {

        for (size_t index = line_offsets.back() + 1; index < block.uncompressed_size; ++index) {

          if (std::isalnum(buffer[index])) {

            ++remain_alpha;

          }

        }

      }

      ExecEnv::log().info("GZBlockDecompression::queueLines; found: {} lines, last line offset: {}, buffer size: {}, remain: {}",
                          line_offsets.size(), line_offsets.empty() ? 0 : line_offsets.back(), block.uncompressed_size, remain_alpha);


    }

  } while (not eof);

  ExecEnv::log().info("Queue Line thread received an EOF notification. Processed: {} blocks", block_count);

}



