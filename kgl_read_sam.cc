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
#include <thread>
#include "kgl_read_sam.h"


namespace kgl = kellerberrin::genome;

// To avoid contention with Python 'log_file', file changed to 'log_file + cpp'
kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : log(sam_read_module_name_, log_file + "cpp")
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

  log.info("Completed processing SAM file");

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

    process_sam_record_.processSamRecord(record_ptr);  // update SNP and Indel structures

    ++counter;

  }

  log.info("Final; Consumer thread processed: {} SAM records", counter);

}