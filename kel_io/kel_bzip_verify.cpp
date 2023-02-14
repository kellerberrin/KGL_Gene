//
// Created by kellerberrin on 21/01/23.
//



#include "kel_bzip.h"

#include "kel_exec_env.h"
#include "kel_utility.h"

#include <fstream>


namespace kel = kellerberrin;


bool kel::BGZReader::verify(const std::string &file_name, bool silent) {

  std::ifstream bgz_file(file_name, std::ios::binary | std::ios::ate);

  if (not bgz_file.good()) {

    ExecEnv::log().warn("BGZStreamIO::verify; Problem opening file: {}", file_name);
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
    BGZHeaderblock header_block;
    bgz_file.read(reinterpret_cast<char*>(&header_block), HEADER_SIZE_);
    file_offset += HEADER_SIZE_;
    // Check the header values.
    if (header_block.block_id_1 != BLOCK_ID1_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, bad header block id 1: {}, expected: {}",
                                            block_count, header_block.block_id_1, BLOCK_ID1_);
      return false;

    }

    if (header_block.block_id_2 != BLOCK_ID2_) {

      if (not silent) ExecEnv::log().error("BGZStreamIO::verify; Block count: {}, bad header block id 2: {}, expected: {}",
                                           block_count, header_block.block_id_2, BLOCK_ID2_);
      return false;

    }

    if (header_block.subfield_id_1 != SUBFIELD_ID1_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                                            block_count, header_block.subfield_id_1, SUBFIELD_ID1_);
      return false;

    }

    if (header_block.subfield_id_2 != SUBFIELD_ID2_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, bad sub-field header block id 1: {}, expected: {}",
                                            block_count, header_block.subfield_id_2, SUBFIELD_ID2_);
      return false;

    }

    if (header_block.length_extra_blocks != EXTRA_LENGTH_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, extra block size: {}, expected: {}",
                                            block_count, header_block.length_extra_blocks, EXTRA_LENGTH_);
      return false;

    }

    // skip the compressed data
    size_t compressed_data_size = header_block.block_size - BLOCK_SIZE_ADJUST_;
    if (compressed_data_size > MAX_UNCOMPRESSED_SIZE_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, Compressed block size: {} exceeds max Ccompressed size: {}",
                                            block_count, compressed_data_size, MAX_UNCOMPRESSED_SIZE_);
      return false;

    }
    bgz_file.seekg(compressed_data_size , std::ios_base::cur);
    file_offset += compressed_data_size;
    total_compressed_size += compressed_data_size;

    // Read the trailer block
    BGZTrailerBlock trailer_block;
    bgz_file.read(reinterpret_cast<char*>(&trailer_block), TRAILER_SIZE_);
    file_offset += TRAILER_SIZE_;
    total_uncompressed_size += trailer_block.uncompressed_size;

    if (trailer_block.uncompressed_size > MAX_UNCOMPRESSED_SIZE_) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, Uncompressed block size: {} exceeds max uncompressed size: {}",
                                            block_count, trailer_block.uncompressed_size, MAX_UNCOMPRESSED_SIZE_);
      return false;

    }

    size_t file_position = static_cast<size_t>(bgz_file.tellg());
    if (file_position != file_offset) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Block count: {}, file tellg: {}, calc file_offset: {}",
                                            block_count, bgz_file.tellg(), file_offset);
      return false;

    }

  } // while.


  // Check EOF string.
  size_t remaining_chars = bgz_file_size - bgz_file.tellg();
  if (remaining_chars != EOF_MARKER_SIZE_) {

    if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; Blocks: Verified {}, EOF Remaining bytes: {}, expected EOF remaining bytes: {}",
                                          block_count, remaining_chars, EOF_MARKER_SIZE_);
    return false;
  }

  uint8_t eof_marker[EOF_MARKER_SIZE_];
  bgz_file.read(reinterpret_cast<char*>(&eof_marker), EOF_MARKER_SIZE_);

  for (size_t index = 0; index < EOF_MARKER_SIZE_; ++index) {

    if (EOF_MARKER_[index] != eof_marker[index]) {

      if (not silent) ExecEnv::log().error( "BGZStreamIO::verify; EOF marker index: {}, EOF marker byte: {}, expected byte: {}",
                                            index, eof_marker[index], EOF_MARKER_[index]);
      return false;

    }

  }

  ExecEnv::log().info("Verified: {}, Blocks: {}, Data Uncompressed: {}, Compressed: {}",
                      file_name,  block_count, total_uncompressed_size, total_compressed_size);

  return true;

}

