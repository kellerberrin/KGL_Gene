// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 18/09/17.
//

#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <regex>
#include <thread>
#include "kgl_consume_sam.h"

namespace kgl = kellerberrin::genome;

// Return false if an unmapped read or problems with cigar.
bool kgl::SAMRecordParser::parseSAM(std::unique_ptr<const std::string>& record_ptr) {

  if (parseSAMFields(record_ptr, false)) {

    if (fastDecodeSAMCigar(record_ptr, sam_fields_[CIGAR_OFFSET])) {

      return true;

    }

  }

  return false;

}

// Efficiently store field offsets and sizes.
bool kgl::SAMRecordParser::parseSAMFields(std::unique_ptr<const std::string>& record_ptr, bool parse_opt_fields)
{

  std::size_t start = 0;
  std::size_t end = record_ptr->find_first_of(SAM_DELIMITER);

  int field_count = 0;

  while (end <= std::string::npos) {

    if (field_count < SAM_FIELD_COUNT) {

      sam_fields_[field_count] = std::pair<std::size_t, std::size_t>(start, end - start);

    } else if (parse_opt_fields) {

      opt_flags_.emplace_back(start, end - start);

    }
    else break;

    ++field_count;

    if (end == std::string::npos)
      break;

    start = end + 1;
    end = record_ptr->find_first_of(SAM_DELIMITER, start);

  }

  if (field_count <  SAM_FIELD_COUNT) {

    log.error("Incorrect field count: {}, SAM record: {}", field_count, *record_ptr);
    return false;

  }

  return (record_ptr->at(sam_fields_[RNAME_OFFSET].first) != UNMAPPED_READ_);

}

bool kgl::SAMRecordParser::decodeSAMCigar( std::unique_ptr<const std::string>& record_ptr
    , const std::pair<std::size_t, std::size_t>& cigar_offset) {


  static std::regex regex("([0-9]+[MIDNSHP=X])");
  static std::sregex_iterator match_end;
  std::smatch match;
  std::string::const_iterator begin = record_ptr->begin() + cigar_offset.first;
  std::string::const_iterator end = begin + cigar_offset.second;

  std::sregex_iterator match_iter(begin, end, regex);

  while (match_iter != match_end) {

    match = *match_iter;
    std::string cigar_item(std::move(match.str()));
    const char cigar_action = cigar_item.back();
    cigar_item.pop_back();
    ContigOffset_t cigar_length = std::stoull(cigar_item);
    std::pair<const char, const ContigOffset_t> cigar(cigar_action, cigar_length);

    cigar_fields_.emplace_back(std::move(cigar));

    match_iter++;

  }

  return not cigar_fields_.empty();

}

bool kgl::SAMRecordParser::fastDecodeSAMCigar( std::unique_ptr<const std::string>& record_ptr
    , const std::pair<std::size_t, std::size_t>& cigar_offset) {

  char cigar_buffer[1000]; // Should big enough
  std::size_t char_index = 0;
  std::size_t buffer_index = 0;

  while (char_index < cigar_offset.second) {

    if (!isdigit(record_ptr->at(cigar_offset.first + char_index))) {

      if (record_ptr->at(cigar_offset.first + char_index) != UNMAPPED_READ_) {

        log.error("Unexpected character in cigar, SAM record: {}", *record_ptr);

      }

      return false;

    }

    while (true) {

      cigar_buffer[buffer_index] = record_ptr->at(cigar_offset.first + char_index);
      if (isalpha(cigar_buffer[buffer_index])) break;
      ++char_index;
      ++buffer_index;
      if (char_index >= cigar_offset.second) {

        return false;

      }

    }
    const char cigar_code = cigar_buffer[buffer_index];
    switch(cigar_code) {

      case 'M':
      case 'X':
      case '=':
      case 'D':
      case 'I':
      case 'S':
      case 'H':
      case 'P':
      case 'N':
        break;

      case UNMAPPED_READ_:
        return false;

      default:
        log.error("Unexpected character in cigar, SAM record: {}", *record_ptr);
        return false;

    }
    cigar_buffer[buffer_index] = '\0';
    ContigOffset_t cigar_length = std::strtoull(cigar_buffer, nullptr, 10);
    cigar_fields_.emplace_back(cigar_code, cigar_length);
    buffer_index = 0;
    ++char_index;

  }

  return true;

}

void kgl::ConsumeMTSAM::finalize() {

  size_t unmapped = unmapped_reads_;
  size_t deleted = delete_nucleotide_;
  size_t inserted = insert_sequence_;
  size_t ignored = other_contig_;
  log.info("Completed processing SAM file; unmapped: {}, deleted: {}, inserted: {}, ignored: {}"
      , unmapped, deleted, inserted, ignored);

}

void kgl::ConsumeMTSAM::consume(std::unique_ptr<const std::string>& record_ptr) {

  SAMRecordParser sam_record_parser(log);

  if (not sam_record_parser.parseSAM(record_ptr)) {

    ++unmapped_reads_;  // Contig is "*"
    return;

  }

  const ContigId_t contig_id(sam_record_parser.getContigId(record_ptr));

  auto matrix_ptr = contig_data_map_.getContig(contig_id);  // Get the contig data block.

  if (matrix_ptr == contig_data_map_.notFound()) {

    ++other_contig_;    // Contig is not registered.
    return;

  }

  ContigArrayMT& contig_block = contig_data_map_.getMatrix(matrix_ptr);

  auto contig_size = contig_block.contigSize();

  ContigOffset_t location = sam_record_parser.getPos(record_ptr);

  if (location >= contig_size) {

    log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}"
        , contig_id, contig_size, location, *record_ptr);
    return;

  }

  ContigOffset_t sam_idx = 0;

  for (auto cigar : sam_record_parser.cigarFields()) {

    switch(cigar.first) {

      case 'M':
      case 'X':
      case '=':
        for (ContigOffset_t cigar_offset = 0; cigar_offset < cigar.second; ++cigar_offset) {

          if (location + cigar_offset >= contig_size) {
            log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {} Cigar Offset: {}; SAM record: {}"
                , contig_id, contig_size, location, cigar_offset, *record_ptr);
            return;
          }
          Nucleotide_t sam_nucleotide = sam_record_parser.getSequenceNucleotide(record_ptr, cigar_offset + sam_idx);
          contig_block.incrementCount(location + cigar_offset, sam_nucleotide);

        }

        sam_idx += cigar.second;
        location += cigar.second;
        break;

      case 'D':

        delete_nucleotide_++;
        for (ContigOffset_t cigar_offset = 0; cigar_offset < cigar.second; ++cigar_offset) {

          contig_block.incrementCount(location + cigar_offset, DELETE_NUCLEOTIDE_);

        }

        location += cigar.second;
        break;

      case 'I':

        insert_sequence_++;
        contig_block.incrementCount(location, INSERT_NUCLEOTIDE_);
        insert_queue_.push(contig_id, location, sam_record_parser.getSubSequence(record_ptr, sam_idx, cigar.second));

        sam_idx += cigar.second;
        break;

      case 'S':
        sam_idx += cigar.second;
        break;

      case 'H':
      case 'P':
        break;

      case 'N':
        location += cigar.second;
        break;

      default:
        log.error("Unknown cigar code {}; SAM record: {}", cigar.first, *record_ptr);
        return;

    }

  }

}

