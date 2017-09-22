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
// Created by kellerberrin on 18/09/17.
//

#ifndef KGL_CONSUME_SAM_H
#define KGL_CONSUME_SAM_H

#include <string>
#include <vector>
#include <queue>
#include "kgl_logging.h"
#include "kgl_mt_data.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


using ConsumerRecordType = ConsumerNumpyRecord;  //Define Local or Numpy
//using ConsumerRecordType = ConsumerLocalRecord;  //Define Local or Numpy
using ConsumerDataMap = ContigDataMap<ConsumerRecordType>;  // A map of contig data blocks.

// Process (consume) SAM records coming from the SAM record reader (the producer).
// This code is multi-threaded, all data structures must be thread safe (see kgl_mt_data.h).

class ConsumeMTSAM {

public:

  explicit ConsumeMTSAM(Logger& logger) : log(logger),  contig_data_map_(logger) {}
  virtual ~ConsumeMTSAM() = default;

  ConsumerDataMap& contigDataMap() { return contig_data_map_; }

  void consume(std::unique_ptr<const std::string>& record_ptr);  // Parse SAM record into fields.
  void finalize();
  void readQuality(unsigned char read_quality) { read_quality_ = read_quality + NUCLEOTIDE_QUALITY_ASCII; }

private:

  Logger& log;                              // Declared First. Emit log messages to console and log file.

  static constexpr unsigned char NUCLEOTIDE_QUALITY_ASCII = 33;  // adjusts nucleotide quality to an ascii value.
  static constexpr unsigned char DEFAULT_MINIMUM_QUALITY = 30 + NUCLEOTIDE_QUALITY_ASCII;   // -10 log10 Pr{ReadError}

  unsigned char read_quality_ = DEFAULT_MINIMUM_QUALITY;
  ConsumerDataMap contig_data_map_;             // Thread safe map of contig data.

  std::atomic<uint64_t> unmapped_reads_{0};     // Contig = "*"
  std::atomic<uint64_t> other_contig_{0};       // SAM records for unregistered contigs (ignored by the code).
  std::atomic<uint64_t> insert_sequence_{0};    // Insertions
  std::atomic<uint64_t> delete_nucleotide_{0};  // Deletions
  std::atomic<uint64_t> rejected_{0};           // Below Quality Reads
  std::atomic<uint64_t> accepted_{0};           // Nucleotides accepted

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_CONSUME_SAM_H
