//
// Created by kellerberrin on 26/02/18.
//

#ifndef KGL_VARIANT_FACTORY_READVCF_IMPL_H
#define KGL_VARIANT_FACTORY_READVCF_IMPL_H


#include <string>
#include <thread>
#include <fstream>
#include "kgl_mt_queue.h"
#include "kgl_exec_env.h"

#include <seqan/vcf_io.h>

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

//template<class ConsumerMT>
//using ConsumerFunctionPtr = void (ConsumerMT::*)(const seqan::VcfRecord& record_ptr);

template<class ConsumerMT>
class VCFReaderMT {

public:

  using ConsumerFunctionPtr = void (ConsumerMT::*)(const seqan::VcfRecord& record_ptr);
  using ConsumerObjPtr = ConsumerMT*;

  explicit VCFReaderMT(ConsumerObjPtr consumer_obj, ConsumerFunctionPtr consumer_fn_ptr) : producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                                                           consumer_obj_ptr_(consumer_obj),
                                                                                           consumer_fn_ptr_(consumer_fn_ptr),
                                                                                           vcf_header_ptr_(std::make_unique<seqan::VcfHeader>())  {
    report_increment_ = 0;

  }
  virtual ~VCFReaderMT() = default;

  void readVCFFile(const std::string &vcf_file_name);

  const seqan::VcfHeader& getHeader() const { return *vcf_header_ptr_; }

  const std::string getContig(int32_t contig_idx) const;

  void setReportIncrement(long report_increment) { report_increment_ = report_increment; }

private:

  BoundedMtQueue<std::unique_ptr<seqan::VcfRecord>> producer_consumer_queue_; // The Producer/Consumer record queue
  ConsumerObjPtr consumer_obj_ptr_;                          // Object to consumer the VCF records.
  ConsumerFunctionPtr consumer_fn_ptr_;                  // Function to consume the VCF records.

  std::shared_ptr<seqan::VcfHeader> vcf_header_ptr_; // the vcf file header record
  std::shared_ptr<seqan::VcfFileIn> vcfIn_ptr_; // Open input file.

  long report_increment_;

  mutable std::mutex mutex_;

  static constexpr long HIGH_TIDE_{100000};          // Maximum BoundedMtQueue size
  static constexpr long LOW_TIDE_{50000};            // Low water mark to begin queueing VCF records

  int consumer_thread_count_{2};                      // Consumer threads (defaults to local CPU cores available)
  static constexpr int MAX_CONSUMER_THREADS_{16};     // Spawning more threads does not increase performance

  // Read the SAM file and queue the record in a BoundedMtQueue.
  void VCFProducer();
  // Call the template SAM consumer class
  void VCFConsumer();

};

template <class ConsumerMT>
const std::string VCFReaderMT<ConsumerMT>::getContig(int32_t contig_idx) const {

  AutoMutex auto_mutex(mutex_); // lock on construction, unlock on destruction

  std::string contig_id = seqan::toCString(contigNames(context(*vcfIn_ptr_))[contig_idx]);

  return contig_id;

}

template <class ConsumerMT>
void VCFReaderMT<ConsumerMT>::readVCFFile( const std::string &vcf_file_name) {

  ExecEnv::log().info("Begin processing VCF file: {}", vcf_file_name);

  vcfIn_ptr_ = std::make_shared<seqan::VcfFileIn>(seqan::toCString(vcf_file_name));

  // Spawn consumer threads.
  consumer_thread_count_ = std::thread::hardware_concurrency();

  // Spawn a maximum of 4 consumers.
  consumer_thread_count_ = consumer_thread_count_ > MAX_CONSUMER_THREADS_? MAX_CONSUMER_THREADS_ : consumer_thread_count_;

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

    long counter = 0;
    while (!seqan::atEnd(*vcfIn_ptr_)) {

      std::unique_ptr<seqan::VcfRecord> record_ptr(std::make_unique<seqan::VcfRecord>());
      seqan::readRecord(*record_ptr, *vcfIn_ptr_);

      producer_consumer_queue_.push(std::move(record_ptr));

      ++counter;

      if (report_increment_ != 0 and counter % report_increment_ == 0) {

        ExecEnv::log().info("Producer thread read: {} VCF records", counter);

      }

    }

    // Enqueue the null eof indicator for each consumer thread.
    for(int i = 0; i < consumer_thread_count_; ++i) {

      std::unique_ptr<seqan::VcfRecord> EOF_INDICATOR{nullptr};  // Enqueued by producer to indicate VCF eof.
      producer_consumer_queue_.push(std::move(EOF_INDICATOR));

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file unexpected I/O exception: {}", e.what());

  }


}

template <class ConsumerMT>
void VCFReaderMT<ConsumerMT>::VCFConsumer() {

  long counter = 0;
  std::unique_ptr<seqan::VcfRecord> record_ptr;
  const std::unique_ptr<seqan::VcfRecord> EOF_INDICATOR{nullptr};

  while (true) {

    producer_consumer_queue_.waitAndPop(record_ptr);

    ++counter;

    if (record_ptr == EOF_INDICATOR) break;  // Eof encountered, terminate processing.

    (consumer_obj_ptr_->*consumer_fn_ptr_)(*record_ptr);

  }

  ExecEnv::log().info("Final; Consumer thread processed: {} VCF records", counter);

}



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
