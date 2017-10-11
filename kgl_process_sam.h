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


// Simple class to bolt the producer, consumer and data together.
using ContigDataBlock = ContigDataMap<ConsumerLocalRecord>;
using LocalConsumer = ConsumeMTSAM<ContigDataBlock>;
using LocalProducer = ProduceMTSAM<LocalConsumer>;

class LocalProcessSam {

public:

  explicit LocalProcessSam(std::shared_ptr<ContigDataBlock> contig_data_ptr, Logger& log) {

    consumer_ptr_ = std::shared_ptr<LocalConsumer>(std::make_shared<LocalConsumer>(log, contig_data_ptr));
    producer_ptr_ = std::unique_ptr<LocalProducer>(std::make_unique<LocalProducer>(log, consumer_ptr_));

  }
  virtual ~LocalProcessSam() = default;

  inline void readSAMFile(const std::string& file_name, unsigned char readQuality = 0) {

    consumer_ptr_->readQuality(readQuality);
    producer_ptr_->readSamFile(file_name);
    consumer_ptr_->finalize();

  }

private:

  std::shared_ptr<LocalConsumer> consumer_ptr_; ;  // Thread safe SAM record consumer
  std::unique_ptr<LocalProducer> producer_ptr_;   //  Thread safe SAM record producer.

};

}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_PROCESS_SAM_H
