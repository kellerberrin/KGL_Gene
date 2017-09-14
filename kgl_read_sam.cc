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

kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : log(sam_read_module_name_, log_file)
                                                                 , producer_consumer_queue_(high_tide_, low_tide_)
                                                                 , process_sam_record_(log)
                                                                 , contig_data_map_(log) {}

void kgl::ProcessSamFile::readSamFile(std::string &file_name) {

  log.info("Begin processing SAM file: {}", file_name);

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency();

  log.info("Spawning: {} Consumer threads to process the SAM file", consumer_thread_count_);

  std::vector<std::thread> consumer_vec;
  for(int i = 0; i < consumer_thread_count_; ++i) {

    consumer_vec.emplace_back(std::thread(&ProcessSamFile::samConsumer,this));

  }

  // Read SAM records and enqueue them.
  samProducer(file_name);

  // Join on the consumer threads

  for(int i = 0; i < consumer_vec.size(); ++i) {

    consumer_vec[i].join();

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

void kgl::ProcessSamFile::parseSAMFields( std::string sam_record
                                        , std::vector<std::string>& sam_fields
                                        , std::vector<std::string>& sam_flags)
{
  constexpr int sam_mandatory_field_count = 11;
  constexpr char sam_delimiter = '\t';
  std::size_t start=0;
  std::size_t end = sam_record.find_first_of(sam_delimiter);

  int field_count = 0;

  while (end <= std::string::npos) {

    ++field_count;

    if (field_count <= sam_mandatory_field_count) {

      sam_fields.emplace_back(sam_record.substr(start, end - start));

    } else {

      sam_flags.emplace_back(sam_record.substr(start, end - start));

    }

    if (end == std::string::npos)
      break;

    start = end + 1;
    end = sam_record.find_first_of(sam_delimiter, start);

  }

  if (sam_fields.size() <  sam_mandatory_field_count) {

    log.error("Incorrect field count: {}, line: {}", sam_fields.size(), sam_record);
    std::exit(EXIT_FAILURE);

  }

}

void kgl::ProcessSamFile::decodeSAMCigar( std::string cigar_string
                                        , std::vector<std::pair<const char, const ContigOffset_t>>& cigar_fields) {

  static std::regex regex("([0-9]+[MIDNSHP=X])");
  static std::sregex_iterator match_end;
  std::smatch match;

  std::sregex_iterator match_iter(cigar_string.begin(), cigar_string.end(), regex);

  while (match_iter != match_end) {

    match = *match_iter;
    std::string cigar_item(std::move(match.str()));
    const char cigar_action = cigar_item.back();
    cigar_item.pop_back();
    ContigOffset_t cigar_length = std::stoull(cigar_item);
    std::pair<const char, const ContigOffset_t> cigar(cigar_action, cigar_length);

    cigar_fields.emplace_back(std::move(cigar));

    match_iter++;

  }

}

void kgl::ProcessSamFile::parseSAMRecord(std::unique_ptr<const std::string>& record_ptr) {

  // Define the SAM record offsets.
  constexpr size_t Qname = 0;
  constexpr size_t Flag = 1;
  constexpr size_t Rname = 2;
  constexpr size_t Pos = 3;
  constexpr size_t Mapquality = 4;
  constexpr size_t Cigar = 5;
  constexpr size_t Rnext = 6;
  constexpr size_t Pnext = 7;
  constexpr size_t Tlen = 8;
  constexpr size_t Sequence = 9;
  constexpr size_t Quality = 10;


  std::vector<std::string> sam_fields;
  std::vector<std::string> sam_flags;
  std::vector<std::pair<const char, const ContigOffset_t>> cigar_fields;

  parseSAMFields(*record_ptr, sam_fields, sam_flags); // Must be called first.
  decodeSAMCigar(sam_fields[Cigar], cigar_fields);  // Must be called first.

  const ContigId_t& contig_id = sam_fields[Rname];  // Alias the contig_id

  if (contig_id == unmapped_read) {

    unmapped_reads_++;
    return;

  }

  ContigMatrixMT& contig_block = contig_data_map_.getContig(contig_id);  // Get the contig data block.

  ContigOffset_t location = std::stoull(sam_fields[Pos]) - 1;    // Subtract 1 for zero offset - SAM usage is offset 1

  auto contig_size = contig_block.contigSize();

  if (location >= contig_size) {

    log.error("Sam record error - Contig: {} sequence size: {} exceeded at position: {}; SAM record: {}"
             , contig_id, contig_size, location, *record_ptr);
    std::exit(EXIT_FAILURE);

  }

  ContigOffset_t sam_idx = 0;

  for (auto cigar : cigar_fields) {

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
          Nucleotide_t sam_nucleotide = sam_fields[Sequence][cigar_offset + sam_idx];
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
          std::string insert_sequence = sam_fields[Sequence].substr(sam_idx, cigar.second);
          insert_sequence_++;
          insert_queue_.push(contig_id, location, insert_sequence);
        }
        sam_idx += cigar.second;
        break;

      case 'S':
        sam_idx += cigar.second;
        break;

      case 'H':
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

