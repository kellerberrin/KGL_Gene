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

#ifndef KGL_PROCESS_SAM_H
#define KGL_PROCESS_SAM_H

#include <memory>
#include "kgl_logging.h"
#include "kgl_read_sam.h"
#include "kgl_consume_sam.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




// Simple class to bolt the Python producer and consumer together.
using LocalConsumer = ConsumeMTSAM<ConsumerLocalRecord>;
using LocalProducer = ProduceMTSAM<LocalConsumer>;

class LocalProcessSam {

public:

  explicit LocalProcessSam(const std::string &log_file, int readQuality) : log(SAM_READ_MODULE_NAME_, log_file) {

    consumer_ptr_ = std::shared_ptr<LocalConsumer>(std::make_shared<LocalConsumer>(log));
    consumer_ptr_->readQuality(readQuality);
    producer_ptr_ = std::unique_ptr<LocalProducer>(std::make_unique<LocalProducer>(log, consumer_ptr_));

  }
  virtual ~LocalProcessSam() = default;

  inline void readSAMFile(const std::string& file_name) {

    producer_ptr_->readSamFile(file_name);
    consumer_ptr_->finalize();

  }

  void insertContig( const ContigId_t& contig_id,
                     const ContigSize_t contig_size) {

    std::unique_ptr<ConsumerLocalRecord> contig_matrix_ptr(std::make_unique<ConsumerLocalRecord>(log,  contig_size));
    contigDataMap().addContigBlock(contig_id, contig_matrix_ptr);

  }


  inline ContigDataMap<ConsumerLocalRecord>& contigDataMap() { return consumer_ptr_->contigDataMap(); }

private:

  Logger log;                              // Must be declared First. Emit log messages to console and log file.
  static constexpr const char* SAM_READ_MODULE_NAME_{"SamRead"};  // Name of this module for the logger

  std::shared_ptr<LocalConsumer> consumer_ptr_; ;  // Thread safe SAM record consumer
  std::unique_ptr<LocalProducer> producer_ptr_;   //  Thread safe SAM record producer.

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_PROCESS_SAM_H
