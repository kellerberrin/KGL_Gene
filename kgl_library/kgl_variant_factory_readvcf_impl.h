//
// Created by kellerberrin on 26/02/18.
//

#ifndef KGL_VARIANT_FACTORY_READVCF_IMPL_H
#define KGL_VARIANT_FACTORY_READVCF_IMPL_H


#include <string>
#include <thread>
#include <fstream>
#include <functional>

#include "kel_lock.h"
#include "kel_mt_queue.h"
#include "kel_exec_env.h"

#include "kgl_variant_file.h"


#include <seqan/vcf_io.h>


namespace kellerberrin::genome {   //  organization::project level namespace



// The vcf record in seqan format and a counter to indicate which vcf record.
template<class ConsumerMT>
class VCFReaderMT {

public:

  using ReaderQueueRecord = std::pair<size_t, std::unique_ptr<VcfRecord>>;
  using ConsumerFunctionPtr = void (ConsumerMT::*)(size_t vcf_record_count, const VcfRecord& record_ptr);
  using ConsumerObjPtr = ConsumerMT*;

  explicit VCFReaderMT(const std::string& vcf_file_name,
                       ConsumerObjPtr consumer_obj,
                       ConsumerFunctionPtr consumer_fn_ptr) : vcf_file_name_(vcf_file_name),
                                                              producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                              consumer_obj_ptr_(consumer_obj),
                                                              consumer_fn_ptr_(consumer_fn_ptr),
                                                              vcf_header_ptr_(std::make_unique<seqan::VcfHeader>())  {

    parseFieldNames(vcf_file_name_);
    vcfIn_ptr_ = std::make_shared<seqan::VcfFileIn>(seqan::toCString(vcf_file_name_));

  }
  virtual ~VCFReaderMT() = default;

  void readVCFFile();

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  [[nodiscard]] const seqan::VcfHeader& readHeader() const;

private:

  const std::string vcf_file_name_;
  BoundedMtQueue<ReaderQueueRecord> producer_consumer_queue_; // The Producer/Consumer record queue
  ConsumerObjPtr consumer_obj_ptr_;                          // Object to consumer the VCF records.
  ConsumerFunctionPtr consumer_fn_ptr_;                  // Function to consume the VCF records.
  VCFParseHeader parseheader_;                          // Get genome and contig information.

  std::shared_ptr<seqan::VcfHeader> vcf_header_ptr_; // the vcf file header record
  std::shared_ptr<seqan::VcfFileIn> vcfIn_ptr_; // Open input file.

  std::unique_ptr<VcfRecord> EOF_INDICATOR_{nullptr};

  mutable std::mutex mutex_;

  int consumer_thread_count_{2};                      // Consumer threads (defaults to local CPU cores available or max)
  static constexpr const int MAX_CONSUMER_THREADS_{6};     // Spawning more threads does not increase performance
  static constexpr const int MIN_CONSUMER_THREADS_{1};     // Need at least 1 consumer thread

  static constexpr const size_t report_increment_{10000};

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records

  [[nodiscard]] ContigId_t getContig(int32_t contig_idx) const;

  // Read the VCF file and queue the record in a BoundedMtQueue.
  void VCFProducer();
  // Call the template VCF consumer class
  void VCFConsumer();

  void parseFieldNames(const std::string& vcf_file_name) { parseheader_.parseHeader(vcf_file_name); }

};



template <class ConsumerMT>
ContigId_t VCFReaderMT<ConsumerMT>::getContig(int32_t contig_idx) const {

  AutoMutex auto_mutex(mutex_); // lock on construction, unlock on destruction

  std::string contig_id = seqan::toCString(contigNames(context(*vcfIn_ptr_))[contig_idx]);

  return contig_id;

}

template <class ConsumerMT>
const seqan::VcfHeader& VCFReaderMT<ConsumerMT>::readHeader() const {

  seqan::readHeader(*vcf_header_ptr_, *vcfIn_ptr_);

  return *vcf_header_ptr_;

}


template <class ConsumerMT>
void VCFReaderMT<ConsumerMT>::readVCFFile() {

  ExecEnv::log().info("Begin processing VCF file: {}", vcf_file_name_);

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency() - 1;

  // Spawn a maximum number of consumers.
  consumer_thread_count_ = consumer_thread_count_ > MAX_CONSUMER_THREADS_? MAX_CONSUMER_THREADS_ : consumer_thread_count_;

  // Spawn a minimum of 1 consumers.
  consumer_thread_count_ = consumer_thread_count_ < MIN_CONSUMER_THREADS_ ? MIN_CONSUMER_THREADS_ : consumer_thread_count_;

  ExecEnv::log().info("Spawning: {} Consumer threads to process the VCF file", consumer_thread_count_);

  std::vector<std::thread> consumer_threads;
  for(int i = 0; i < consumer_thread_count_; ++i) {

    consumer_threads.emplace_back(&VCFReaderMT::VCFConsumer, this);

  }

  // Read VCF records and enqueue them.
  VCFProducer();

  // Join on the consumer threads
  for(auto& thread : consumer_threads) {

    thread.join();

  }

}

// Read the SAM file and queue the records.
template<class ConsumerMT>
void VCFReaderMT<ConsumerMT>::VCFProducer() {

  try {

    seqan::readHeader(*vcf_header_ptr_, *vcfIn_ptr_);

    size_t counter = 1;
    while (!seqan::atEnd(*vcfIn_ptr_)) {

      seqan::VcfRecord vcf_record;
      seqan::readRecord(vcf_record, *vcfIn_ptr_);

      std::unique_ptr<VcfRecord> record_ptr(std::make_unique<VcfRecord>(std::move(vcf_record), getContig(vcf_record.rID)));

      ReaderQueueRecord queue_record(counter, std::move(record_ptr));

      producer_consumer_queue_.push(std::move(queue_record));

      if (counter % report_increment_ == 0) {

        ExecEnv::log().info("Producer thread read: {} VCF records", counter);

      }

      ++counter;

    }

    // Enqueue the null eof indicator for each consumer thread.
    for(int i = 0; i < consumer_thread_count_; ++i) {

      ReaderQueueRecord queue_record(counter, std::move(EOF_INDICATOR_));

      producer_consumer_queue_.push(std::move(queue_record));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected I/O exception: {}", e.what());

  }


}

template <class ConsumerMT>
void VCFReaderMT<ConsumerMT>::VCFConsumer() {

  long counter = 0;
  ReaderQueueRecord queue_record;

  // Loop until EOF.
  while (true) {

    // Dequeue the vcf record.
    producer_consumer_queue_.waitAndPop(queue_record);

    ++counter;

    // Terminate on EOF
    if (queue_record.second == EOF_INDICATOR_) break;  // Eof encountered, terminate processing.

    // Call the consumer object with the dequeued record.
    std::invoke(consumer_fn_ptr_, consumer_obj_ptr_, queue_record.first, *queue_record.second);

  }

  ExecEnv::log().info("Final; Consumer thread processed: {} VCF records", counter);

}


}   // end namespace


#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
