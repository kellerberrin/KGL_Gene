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
#include <fstream>
#include <iostream>
#include <seqan/bam_io.h>
#include "kgl_read_sam.h"


namespace kgl = kellerberrin::genome;

// To avoid contention with Python 'log_file', 'libread_sam' logs to 'log_file + cpp'
kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : log("ReadSam", log_file + "cpp") {}

void kgl::ProcessSamFile::readSamFile(std::string &file_name) {

  log.info("'libread_sam': Begin processing SAM file");

  // Open input file, BamFileIn can read SAM and BAM files.
  seqan::BamFileIn bamFileIn;

  if (!open(bamFileIn, seqan::toCString(file_name))) {
    log.error("I/O error; could not open file: {}", file_name);
    return;
  }

  try {
    // Copy header.
    seqan::BamHeader header;
    readHeader(header, bamFileIn);
    long counter = 0;

    // Copy records.
    seqan::BamAlignmentRecord record;
    while (not atEnd(bamFileIn)) {

      readRecord(record, bamFileIn);
      ++counter;

    }

    log.info("Processed: {} records", counter);

  }
  catch (seqan::Exception const &e) {
    log.error("seqan library exception: {}", e.what());
    return;
  }

  log.info("'libread_sam': Completed processing SAM file");

  return;
}
