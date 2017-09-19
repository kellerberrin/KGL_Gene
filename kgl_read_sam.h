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
// Created by kellerberrin on 8/09/17.
//

#ifndef KGL_READ_SAM_H
#define KGL_READ_SAM_H

#include <string>
#include <thread>
#include "kgl_logging.h"
#include "kgl_mt_queue.h"
#include "kgl_consume_sam.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class ProduceMTSAM {

public:

  explicit ProduceMTSAM(Logger& logger, ConsumeMTSAM& consumer) : log(logger)
                                                                , producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_)
                                                                , consumer_(consumer) {}
  virtual ~ProduceMTSAM() = default;

  void readSamFile(const std::string &file_name);

private:

  Logger& log;                              // Declared First. Emit log messages to console and log file.
  BoundedMtQueue<std::unique_ptr<const std::string>> producer_consumer_queue_; // The Producer/Consumer record queue
  ConsumeMTSAM& consumer_;                  // Consume the SAM records.

  static constexpr const char* EOF_INDICATOR_{"<<EOF>>"};  // Enqueued by producer to indicate SAM eof.
  static constexpr long REPORT_INCREMENT_{500000};    // Frequency to emit SAM progress messages
  static constexpr long HIGH_TIDE_{1000000};          // Maximum BoundedMtQueue size
  static constexpr long LOW_TIDE_{500000};            // Low water mark to begin queueing SAM records

  int consumer_thread_count_{4};                      // Consumer threads (defaults to local CPU cores available)

  // Read the SAM file and queue the record in a BoundedMtQueue.
  void samProducer(const std::string &file_name);
  // Call the template SAM consumer class
  void samConsumer();


};

}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_READ_SAM_H
