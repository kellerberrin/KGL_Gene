//
// Created by kellerberrin on 16/4/20.
//

#include "kgl_variant_factory_readvcf_impl.h"


namespace kgl = kellerberrin::genome;


void kgl::VCFReaderMT::readHeader() {

  processVCFHeader(vcf_io_.VCFReadHeader());

}


void kgl::VCFReaderMT::readVCFFile() {


  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency() - 1;

  // Spawn a maximum number of consumers.
  consumer_thread_count_ = consumer_thread_count_ > MAX_CONSUMER_THREADS_? MAX_CONSUMER_THREADS_ : consumer_thread_count_;

  // Spawn a minimum of 1 consumers.
  consumer_thread_count_ = consumer_thread_count_ < MIN_CONSUMER_THREADS_ ? MIN_CONSUMER_THREADS_ : consumer_thread_count_;

  ExecEnv::log().info("Parse Header VCF file: {}", vcf_io_.getFileName());

  readHeader();

  ExecEnv::log().info("Begin processing VCF file: {}", vcf_io_.getFileName());

  ExecEnv::log().info("Spawning: {} Consumer threads to process the VCF file", consumer_thread_count_);

  std::vector<std::thread> consumer_threads;
  for(size_t i = 0; i < consumer_thread_count_; ++i) {

    consumer_threads.emplace_back(&VCFReaderMT::VCFConsumer, this);

  }

  // start reading records asynchronously.
  vcf_io_.commenceVCFIO(consumer_thread_count_);

  // Join on the consumer threads
  for(auto& thread : consumer_threads) {

    thread.join();

  }

}


void kgl::VCFReaderMT::VCFConsumer() {

  // Loop until EOF.
  size_t final_count = 0;
  while (true) {

    // Dequeue the vcf record.
    QueuedVCFRecord vcf_record = vcf_io_.readVCFRecord();

    // Terminate on EOF
    if (not vcf_record) break;  // Eof encountered, terminate processing.

    if (vcf_record.value().second->contig_id.empty()) {

      ExecEnv::log().error("Empty VCF_record encountered; consumer thread terminates.");
      break;

    }

    // Call the consumer object with the dequeued record.
    ProcessVCFRecord(vcf_record.value().first, *vcf_record.value().second);
    final_count = vcf_record.value().first;

  }

  ExecEnv::log().info("Final; Consumer thread processed: {} VCF records", final_count);

}


