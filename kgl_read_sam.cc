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
//
//
// Created by kellerberrin on 26/08/17.
//

#include <string>
#include <memory>
#include <fstream>
#include <iostream>
#include <regex>
#include <thread>
#include "kgl_read_sam.h"

namespace kgl = kellerberrin::genome;


bool kgl::SAMRecordParser::Parse(std::unique_ptr<const std::string>& record_ptr) {

  parseSAMFields(record_ptr, false);
  if (record_ptr->at(sam_fields_[RNAME_OFFSET].first) != unmapped_read) {

    if (fastDecodeSAMCigar(record_ptr, sam_fields_[CIGAR_OFFSET])) {

      return true;

    }

  }

  return false;

}

// Efficiently store field offsets and sizes.
void kgl::SAMRecordParser::parseSAMFields(std::unique_ptr<const std::string>& record_ptr, bool parse_opt_fields)
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
    std::exit(EXIT_FAILURE);

  }

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

  return true;

}

bool kgl::SAMRecordParser::fastDecodeSAMCigar( std::unique_ptr<const std::string>& record_ptr
                                             , const std::pair<std::size_t, std::size_t>& cigar_offset) {

  char cigar_buffer[1000];
  std::size_t char_index = 0;
  std::size_t buffer_index = 0;

  while (char_index < cigar_offset.second) {

    if (!isdigit(record_ptr->at(cigar_offset.first + char_index))) {

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

      default:
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

kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : log(sam_read_module_name_, log_file)
                                                                 , producer_consumer_queue_(high_tide_, low_tide_)
                                                                 , contig_data_map_(log) {}

void kgl::ProcessSamFile::readSamFile(std::string &file_name) {

  log.info("Begin processing SAM file: {}", file_name);

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency();

  log.info("Spawning: {} Consumer threads to process the SAM file", consumer_thread_count_);

  std::vector<std::thread> consumer_threads;
  for(int i = 0; i < consumer_thread_count_; ++i) {

    consumer_threads.emplace_back(std::thread(&ProcessSamFile::samConsumer,this));

  }

  // Read SAM records and enqueue them.
  samProducer(file_name);

  // Join on the consumer threads
  for(auto& consumer : consumer_threads) {

    consumer.join();

  }

  size_t unmapped = unmapped_reads_;
  size_t deleted = delete_nucleotide_;
  size_t inserted = insert_sequence_;
  log.info("Completed processing SAM file; unmapped: {}, deleted: {}, inserted: {}", unmapped, deleted, inserted);

}

// Read the SAM file and queue the records.
void kgl::ProcessSamFile::samProducer(std::string &file_name) {

  std::ifstream sam_file;

  // Open input file.

  sam_file.open(file_name);

  if (not sam_file.good()) {
    log.error("I/O error; could not open SAM file: {}, program exits.", file_name);
    std::exit(EXIT_FAILURE);
  }

  try {

    long counter = 0;

    while (true) {

      std::unique_ptr<std::string> record_ptr(std::make_unique<std::string>());

      if (std::getline(sam_file, *record_ptr).eof()) break;

      if ((*record_ptr)[0] == '@') continue;   // ignore header records.

      producer_consumer_queue_.push(std::move(record_ptr));

      ++counter;

      if (counter % report_increment_ == 0) {

        log.info("Producer thread read: {} SAM records", counter);

      }

    }

    // Enqueue an eof indicator for each consumer thread.
    for(int i = 0; i < consumer_thread_count_; ++i) {

      std::unique_ptr<std::string> eof_record_ptr(new std::string(eof_indicator_));
      producer_consumer_queue_.push(std::move(eof_record_ptr));

    }

    sam_file.close();

    log.info("Final; Producer thread read: {} SAM records", counter);

  }
  catch (std::exception const &e) {
    log.error("SAM file: {}, unexpected I/O exception: {}, program exits.", file_name, e.what());
    std::exit(EXIT_FAILURE);
  }

}

// Multiple threads; dequeue from the BoundedMtQueue and process.
void kgl::ProcessSamFile::samConsumer() {

  long counter = 0;
  std::unique_ptr<const std::string> record_ptr;
  const std::string eof_record(eof_indicator_);

  while (true) {

    producer_consumer_queue_.waitAndPop(record_ptr);

    if (*record_ptr == eof_record) break;  // Eof encountered, terminate processing.

    parseSAMRecord(record_ptr); // update SNP and Indel structures

    ++counter;

  }

  log.info("Final; Consumer thread processed: {} SAM records", counter);

}

void kgl::ProcessSamFile::parseSAMRecord(std::unique_ptr<const std::string>& record_ptr) {

  SAMRecordParser sam_record_parser(log);

  if (not sam_record_parser.Parse(record_ptr)) {

    unmapped_reads_++;
    return;

  }

  const ContigId_t contig_id(sam_record_parser.getContigId(record_ptr));

  ContigMatrixMT& contig_block = contig_data_map_.getContig(contig_id);  // Get the contig data block.

  ContigOffset_t location = sam_record_parser.getPos(record_ptr);

  auto contig_size = contig_block.contigSize();

  if (location >= contig_size) {

    log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}"
             , contig_id, contig_size, location, *record_ptr);
    std::exit(EXIT_FAILURE);

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
            std::exit(EXIT_FAILURE);
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

          contig_block.incrementCount(location + cigar_offset, delete_nucleotide);

        }

        location += cigar.second;
        break;

      case 'I':

        if (location >= contig_size) {
          log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}"
              , contig_id, contig_size, location, *record_ptr);
          std::exit(EXIT_FAILURE);
        }
        contig_block.incrementCount(location, insert_nucleotide);
        { // Program block to scope insert_sequence variable

          insert_sequence_++;
          insert_queue_.push(contig_id, location, sam_record_parser.getSubSequence(record_ptr, sam_idx, cigar.second));
        }
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
        std::exit(EXIT_FAILURE);

    }

  }

}

