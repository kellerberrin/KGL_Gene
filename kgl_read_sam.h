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

template<class SAMConsumerMT>
class ProduceMTSAM {

public:

  explicit ProduceMTSAM(Logger& logger, std::shared_ptr<SAMConsumerMT>& consumer_ptr)
      : log(logger), producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_) {

    consumer_ptr_ = consumer_ptr;

  }
  virtual ~ProduceMTSAM() = default;

  void readSamFile(const std::string &file_name);

private:

  Logger& log;                              // Declared First. Emit log messages to console and log file.
  BoundedMtQueue<std::unique_ptr<const std::string>> producer_consumer_queue_; // The Producer/Consumer record queue
  std::shared_ptr<SAMConsumerMT> consumer_ptr_;                  // Consume the SAM records.

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

template <class SAMConsumerMT>
void ProduceMTSAM<SAMConsumerMT>::readSamFile( const std::string &file_name) {

  log.info("Begin processing SAM file: {}", file_name);

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency();

  log.info("Spawning: {} Consumer threads to process the SAM file", consumer_thread_count_);

  std::vector<std::thread> consumer_threads;
  for(int i = 0; i < consumer_thread_count_; ++i) {

    consumer_threads.emplace_back(&ProduceMTSAM::samConsumer, this);

  }

  // Read SAM records and enqueue them.
  samProducer(file_name);

  // Join on the consumer threads
  for(auto& thread : consumer_threads) {

    thread.join();

  }

}

// Read the SAM file and queue the records.
template<class SAMConsumerMT>
void ProduceMTSAM<SAMConsumerMT>::samProducer(const std::string &file_name) {

  std::ifstream sam_file;

  // Open input file.

  sam_file.open(file_name);

  if (not sam_file.good()) {

    log.critical("I/O error; could not open SAM file: {}", file_name);

  }

  try {

    long counter = 0;

    while (true) {

      std::unique_ptr<std::string> record_ptr(std::make_unique<std::string>());

      if (std::getline(sam_file, *record_ptr).eof()) break;

      if ((*record_ptr)[0] == '@') continue;   // ignore header records.

      producer_consumer_queue_.push(std::move(record_ptr));

      ++counter;

      if (counter % REPORT_INCREMENT_ == 0) {

        log.info("Producer thread read: {} SAM records", counter);

      }

    }

    // Enqueue an eof indicator for each consumer thread.
    for(int i = 0; i < consumer_thread_count_; ++i) {

      std::unique_ptr<std::string> eof_record_ptr(new std::string(EOF_INDICATOR_));
      producer_consumer_queue_.push(std::move(eof_record_ptr));

    }

    sam_file.close();

    log.info("Final; Producer thread read: {} SAM records", counter);

  }
  catch (std::exception const &e) {

    log.critical("SAM file: {}, unexpected I/O exception: {}", file_name, e.what());

  }

}

template <class SAMConsumerMT>
void ProduceMTSAM<SAMConsumerMT>::samConsumer() {

  long counter = 0;
  std::unique_ptr<const std::string> record_ptr;
  const std::string eof_record(EOF_INDICATOR_);

  while (true) {

    producer_consumer_queue_.waitAndPop(record_ptr);

    if (*record_ptr == eof_record) break;  // Eof encountered, terminate processing.

    consumer_ptr_->consume(record_ptr);

    ++counter;

  }

  log.info("Final; Consumer thread processed: {} SAM records", counter);

}

}   // namespace genome
}   // namespace kellerberrin

#endif // KGL_READ_SAM_H
