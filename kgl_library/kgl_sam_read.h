//
// Created by kellerberrin on 31/10/17.
//

#ifndef KGL_SAM_READ_H
#define KGL_SAM_READ_H


#include <string>
#include <thread>
#include <fstream>
#include "kgl_mt_queue.h"
#include "kgl_exec_env.h"

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

template<class SAMConsumerMT>
class ProduceMTSAM {

public:

  explicit ProduceMTSAM(std::shared_ptr<SAMConsumerMT>& consumer_ptr) : producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_) {

    consumer_ptr_ = consumer_ptr;

  }
  virtual ~ProduceMTSAM() = default;

  void readSamFile(const std::string &file_name);

private:

  BoundedMtQueue<std::unique_ptr<const std::string>> producer_consumer_queue_; // The Producer/Consumer record queue
  std::shared_ptr<SAMConsumerMT> consumer_ptr_;                  // Consume the SAM records.

  static constexpr const char* EOF_INDICATOR_{"<<EOF>>"};  // Enqueued by producer to indicate SAM eof.
  static constexpr long REPORT_INCREMENT_{1000000};    // Frequency to emit SAM progress messages
  static constexpr long HIGH_TIDE_{1000000};          // Maximum BoundedMtQueue size
  static constexpr long LOW_TIDE_{500000};            // Low water mark to begin queueing SAM records

  int consumer_thread_count_{2};                      // Consumer threads (defaults to local CPU cores available)
  static constexpr int MAX_CONSUMER_THREADS_{4};     // Spawning more threads does not increase performance

  // Read the SAM file and queue the record in a BoundedMtQueue.
  void samProducer(const std::string &file_name);
  // Call the template SAM consumer class
  void samConsumer();


};

template <class SAMConsumerMT>
void ProduceMTSAM<SAMConsumerMT>::readSamFile( const std::string &file_name) {

  ExecEnv::log().info("Begin processing SAM file: {}", file_name);

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency();

  // Spawn a maximum of 4 consumers.
  consumer_thread_count_ = consumer_thread_count_ > MAX_CONSUMER_THREADS_? MAX_CONSUMER_THREADS_ : consumer_thread_count_;

  ExecEnv::log().info("Spawning: {} Consumer threads to process the SAM file", consumer_thread_count_);

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

    ExecEnv::log().critical("I/O error; could not open SAM file: {}", file_name);

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

        ExecEnv::log().info("Producer thread read: {} SAM records", counter);

      }

    }

    // Enqueue an eof indicator for each consumer thread.
    for(int i = 0; i < consumer_thread_count_; ++i) {

      std::unique_ptr<std::string> eof_record_ptr(new std::string(EOF_INDICATOR_));
      producer_consumer_queue_.push(std::move(eof_record_ptr));

    }

    sam_file.close();

    ExecEnv::log().info("Final; Producer thread read: {} SAM records", counter);

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("SAM file: {}, unexpected I/O exception: {}", file_name, e.what());

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

  ExecEnv::log().info("Final; Consumer thread processed: {} SAM records", counter);

}

}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_SAM_READ_H
