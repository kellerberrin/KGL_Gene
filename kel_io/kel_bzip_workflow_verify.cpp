//
// Created by kellerberrin on 17/01/23.
//


#include "kel_bzip_workflow.h"
#include "kel_exec_env.h"

#include <fstream>


namespace kel = kellerberrin;


bool kel::BGZStreamIO::verify(const std::string &file_name, bool silent) {

  std::ifstream bgz_file(file_name, std::ios::binary | std::ios::ate);

  if (not bgz_file.good()) {

    ExecEnv::log().warn("Problem opening file: {}", file_name);
    return false;

  }

  try {

    // Read file size.
    const size_t bgz_file_size = static_cast<size_t>(bgz_file.tellg());
    ExecEnv::log().info("Verifying bgz file: {}:, Size: {}", file_name, bgz_file_size);
    // Reset the file stream pointer to the file beginning.
    bgz_file.seekg(0);

    size_t file_offset{0};
    size_t block_count{0};
    size_t total_compressed_size{0};
    size_t total_uncompressed_size{0};
    if (bgz_file_size <= EOF_MARKER_SIZE_) {

      if (not silent) ExecEnv::log().error("File: {} too small, size: {}", file_name, bgz_file_size);
      return false;

    }
    while (file_offset < (bgz_file_size - EOF_MARKER_SIZE_) and not bgz_file.eof()) {

      if (not bgz_file.good()) {

        if (not silent) ExecEnv::log().error("Stream error while verifying file: {}", file_name);
        return false;

      }

      ++block_count;
      // Read header.
      BGZHeaderblock header_block;
      bgz_file.read(reinterpret_cast<char*>(&header_block), HEADER_SIZE_);
      if (not bgz_file.good()) {

        if (not silent) ExecEnv::log().warn("Error reading file: {}", file_name);
        return false;

      }

      file_offset += HEADER_SIZE_;
      // Check the header values.
      if (header_block.block_id_1 != BLOCK_ID1_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, bad header block id 1: {}, expected: {}",
                                              block_count, header_block.block_id_1, BLOCK_ID1_);
        return false;

      }

      if (header_block.block_id_2 != BLOCK_ID2_) {

        if (not silent) ExecEnv::log().error("Block count: {}, bad header block id 2: {}, expected: {}",
                                             block_count, header_block.block_id_2, BLOCK_ID2_);
        return false;

      }

      if (header_block.subfield_id_1 != SUBFIELD_ID1_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                                              block_count, header_block.subfield_id_1, SUBFIELD_ID1_);
        return false;

      }

      if (header_block.subfield_id_2 != SUBFIELD_ID2_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, bad sub-field header block id 2: {}, expected: {}",
                                              block_count, header_block.subfield_id_2, SUBFIELD_ID2_);
        return false;

      }

      if (header_block.length_extra_blocks != EXTRA_LENGTH_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, extra block size: {}, expected: {}",
                                              block_count, header_block.length_extra_blocks, EXTRA_LENGTH_);
        return false;

      }

      // skip the compressed data
      if (header_block.block_size <= BLOCK_SIZE_ADJUST_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, invalid block size: {}, adjust: {}",
                                              block_count, header_block.block_size, BLOCK_SIZE_ADJUST_);
        return false;

      }
      size_t compressed_data_size = header_block.block_size - BLOCK_SIZE_ADJUST_;
      if (compressed_data_size > MAX_UNCOMPRESSED_SIZE_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, Compressed block size: {} exceeds max compressed size: {}",
                                              block_count, compressed_data_size, MAX_UNCOMPRESSED_SIZE_);
        return false;

      }
      bgz_file.seekg(compressed_data_size , std::ios_base::cur);
      file_offset += compressed_data_size;
      total_compressed_size += compressed_data_size;

      // Read the trailer block
      if (not bgz_file.good()) {

        if (not silent) ExecEnv::log().error("Block count: {}, stream error before trailer read", block_count);
        return false;

      }
      BGZTrailerBlock trailer_block;
      bgz_file.read(reinterpret_cast<char*>(&trailer_block), TRAILER_SIZE_);
      if (not bgz_file.good()) {

        if (not silent) ExecEnv::log().error("Block count: {}, error reading trailer block", block_count);
        return false;

      }
      file_offset += TRAILER_SIZE_;
      total_uncompressed_size += trailer_block.uncompressed_size;

      if (trailer_block.uncompressed_size > MAX_UNCOMPRESSED_SIZE_) {

        if (not silent) ExecEnv::log().error( "Block count: {}, Uncompressed block size: {} exceeds max uncompressed size: {}",
                                              block_count, trailer_block.uncompressed_size, MAX_UNCOMPRESSED_SIZE_);
        return false;

      }

      size_t file_position = static_cast<size_t>(bgz_file.tellg());
      if (file_position != file_offset) {

        if (not silent) {

          ExecEnv::log().error( "Block count: {}, file tellg: {}, calc file_offset: {}",
                                block_count, file_position, file_offset);

        }
        return false;

      }

    } // end while.


    // Check stream state before using tellg().
    if (not bgz_file.good()) {

      if (not silent) ExecEnv::log().error("Stream error before EOF marker check, file: {}", file_name);
      return false;

    }
    // Check EOF string.
    size_t remaining_chars = bgz_file_size - static_cast<size_t>(bgz_file.tellg());
    if (remaining_chars != EOF_MARKER_SIZE_) {

      if (not silent) ExecEnv::log().error( "Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                                            block_count, remaining_chars, EOF_MARKER_SIZE_);
      return false;
    }

    uint8_t eof_marker[EOF_MARKER_SIZE_];
    bgz_file.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);
    if (not bgz_file.good()) {

      if (not silent) ExecEnv::log().error("Error reading EOF marker from file: {}", file_name);
      return false;

    }

    for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

      if (EOF_MARKER_[index] != eof_marker[index]) {

        if (not silent) ExecEnv::log().error( "EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                                              index, eof_marker[index], EOF_MARKER_[index]);
        return false;

      }

    }

    ExecEnv::log().info("Verified: {}, Blocks: {}, Data Uncompressed: {}, Compressed: {}",
                        file_name,  block_count, total_uncompressed_size, total_compressed_size);

    return true;

  }
  catch (std::exception const &e) {

    if (not silent) ExecEnv::log().error("BGZStreamIO::verify; exception verifying file: {} error: {}", file_name, e.what());
    return false;

  }

}
