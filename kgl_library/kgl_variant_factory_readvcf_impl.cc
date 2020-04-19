//
// Created by kellerberrin on 16/4/20.
//

#include "kgl_variant_factory_readvcf_impl.h"


namespace kgl = kellerberrin::genome;


void kgl::VCFReaderMT::readHeader() {

  processVCFHeader(vcf_io_ptr_->VCFReadHeader());

}


void kgl::VCFReaderMT::readVCFFile() {


  ExecEnv::log().info("Parse Header VCF file: {}", vcf_io_ptr_->fileName());

  readHeader();

  ExecEnv::log().info("Begin processing VCF file: {}", vcf_io_ptr_->fileName());

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency() - 1;

  // Spawn a maximum number of consumers.
  consumer_thread_count_ = consumer_thread_count_ > MAX_CONSUMER_THREADS_? MAX_CONSUMER_THREADS_ : consumer_thread_count_;

  // Spawn a minimum of 1 consumers.
  consumer_thread_count_ = consumer_thread_count_ < MIN_CONSUMER_THREADS_ ? MIN_CONSUMER_THREADS_ : consumer_thread_count_;

  ExecEnv::log().info("Spawning: {} Consumer threads to process the VCF file", consumer_thread_count_);

  std::vector<std::thread> consumer_threads;
  for(size_t i = 0; i < consumer_thread_count_; ++i) {

    consumer_threads.emplace_back(&VCFReaderMT::VCFConsumer, this);

  }

  // start reading records asynchronously.

  vcf_io_ptr_->commenceIO();

  // Read VCF records and enqueue them.
  VCFProducer();

  // Join on the consumer threads
  for(auto& thread : consumer_threads) {

    thread.join();

  }

}

// Read the VCF file and queue the records.
void kgl::VCFReaderMT::VCFProducer() {

  try {

    size_t counter = 1;
    while (true) {

      std::unique_ptr<VcfRecord> record_ptr = vcf_io_ptr_->readVCFRecord();

      if (not record_ptr) break;

      ReaderQueueRecord queue_record(counter, std::move(record_ptr));

      producer_consumer_queue_.push(std::move(queue_record));

      if (counter % report_increment_ == 0) {

        ExecEnv::log().info("Producer thread read: {} VCF records", counter);

      }

      ++counter;

    }

    // Enqueue the null eof indicator for each consumer thread.
    for(size_t i = 0; i < consumer_thread_count_; ++i) {

      ReaderQueueRecord queue_record(counter, VcfRecord::EOF_RECORD());

      producer_consumer_queue_.push(std::move(queue_record));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected I/O exception: {}", e.what());

  }


}

void kgl::VCFReaderMT::VCFConsumer() {

  long counter = 0;
  ReaderQueueRecord queue_record;

  // Loop until EOF.
  while (true) {

    // Dequeue the vcf record.
    producer_consumer_queue_.waitAndPop(queue_record);

    ++counter;

    // Terminate on EOF
    if (queue_record.second == VcfRecord::EOF_RECORD()) break;  // Eof encountered, terminate processing.

    // Call the consumer object with the dequeued record.
    ProcessVCFRecord(queue_record.first, *queue_record.second);

  }

  ExecEnv::log().info("Final; Consumer thread processed: {} VCF records", counter);

}


