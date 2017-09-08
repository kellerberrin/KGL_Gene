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
#include <seqan/bam_io.h>
#include "kgl_read_sam.h"
#include "kgl_mt_queue.h"


namespace kgl = kellerberrin::genome;

// To avoid contention with Python 'log_file', 'libread_sam' logs to 'log_file + cpp'
kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : log(sam_read_module_name_, log_file + "cpp") {}


void kgl::ProcessSamFile::readSamFile(std::string &file_name) {

  std::ifstream sam_file;

  log.info("Begin processing SAM file: {}", file_name);

  // Open input file, BamFileIn can read SAM and BAM filets.

  sam_file.open(file_name);

  if (not sam_file.good()) {
    log.error("I/O error; could not open file: {}", file_name);
    return;
  }

  try {

    // The bounded multithreaded queue has a maximum of 1000000 elements and a low tide
    // of 500000 when the producer(s) can start pushing elements after a high tide event.
    // This stops excessive memory usage (and swapping) if the producer(s) can queue records
    // faster than consumer(s) can remove them.
    kgl::BoundedMtQueue<std::unique_ptr<std::string>> record_vec(1000000, 500000);
    long counter = 0;

    while (true) {

      std::unique_ptr<std::string> record_ptr(std::make_unique<std::string>());

      if (std::getline(sam_file, *record_ptr).eof()) break;

      if ((*record_ptr)[0] == '@') continue;   // ignore header records.

      record_vec.push(std::move(record_ptr));

      ++counter;

      if (counter % report_increment_ == 0) {

        log.info("Processed: {} records", counter);

      }

    }

    sam_file.close();

    log.info("Final; Processed: {} records", counter);

  }
  catch (std::exception const &e) {
    log.error("SAM file: {} io exception: {}", file_name, e.what());
    return;
  }

}
