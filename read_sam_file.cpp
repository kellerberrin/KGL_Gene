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
#include "read_sam_file.h"


namespace kgl = kellerberrin::genome;

kgl::ProcessSamFile::ProcessSamFile(const std::string& log_file) : logger_(log_file ) {}

void kgl::ProcessSamFile::zreadSamFile(std::string &file_name) {

  long counter = 0;
  std::string samrecord;
  std::ifstream samfile;

  samfile.open(file_name);

  while (not std::getline(samfile, samrecord).eof()) {


    counter++;

  }

  samfile.close();

  std::cout << "SAM records:" << counter << std::endl;

}

using namespace seqan;

void kgl::ProcessSamFile::readSamFile(std::string &file_name) {

  // Open input file, BamFileIn can read SAM and BAM files.
  BamFileIn bamFileIn;

  if (!open(bamFileIn, toCString(file_name))) {
    std::cerr << "ERROR: Could not open " << file_name << std::endl;
    return;
  }

  try {
    // Copy header.
    BamHeader header;
    readHeader(header, bamFileIn);
    long counter = 0;

    // Copy records.
    BamAlignmentRecord record;
    while (not atEnd(bamFileIn)) {

      readRecord(record, bamFileIn);
      ++counter;

    }

    std::cout << "Seqan Read SAM Records:" << counter << std::endl;

  }
  catch (Exception const &e) {
    std::cout << "ERROR: " << e.what() << std::endl;
    return;
  }


  return;
}
