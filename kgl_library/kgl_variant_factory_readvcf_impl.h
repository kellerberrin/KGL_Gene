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

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace bt = boost;

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

//template<class ConsumerMT>
//using ConsumerFunctionPtr = void (ConsumerMT::*)(const seqan::VcfRecord& record_ptr);

template<class ConsumerMT>
class VCFReaderMT {

public:

  using ConsumerFunctionPtr = void (ConsumerMT::*)(const seqan::VcfRecord& record_ptr);
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
    report_increment_ = 0;

  }
  virtual ~VCFReaderMT() = default;

  const seqan::VcfHeader& readHeader() const;

  void readVCFFile();

  const std::vector<std::string>& getFieldNames() const { return field_names_; }

  const ContigId_t getContig(int32_t contig_idx) const;

  void setReportIncrement(long report_increment) { report_increment_ = report_increment; }

private:

  const std::string vcf_file_name_;
  BoundedMtQueue<std::unique_ptr<seqan::VcfRecord>> producer_consumer_queue_; // The Producer/Consumer record queue
  ConsumerObjPtr consumer_obj_ptr_;                          // Object to consumer the VCF records.
  ConsumerFunctionPtr consumer_fn_ptr_;                  // Function to consume the VCF records.
  std::vector<std::string> field_names_;                // Field (genome) names for the VCF record

  std::shared_ptr<seqan::VcfHeader> vcf_header_ptr_; // the vcf file header record
  std::shared_ptr<seqan::VcfFileIn> vcfIn_ptr_; // Open input file.

  long report_increment_;

  mutable std::mutex mutex_;

  static constexpr long HIGH_TIDE_{100000};          // Maximum BoundedMtQueue size
  static constexpr long LOW_TIDE_{50000};            // Low water mark to begin queueing VCF records

  int consumer_thread_count_{2};                      // Consumer threads (defaults to local CPU cores available)
  static constexpr int MAX_CONSUMER_THREADS_{16};     // Spawning more threads does not increase performance
  static constexpr int MIN_CONSUMER_THREADS_{1};     // Need at least 1 producer thread

  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.

  // Read the VCF file and queue the record in a BoundedMtQueue.
  void VCFProducer();
  // Call the template VCF consumer class
  void VCFConsumer();

  void parseFieldNames(const std::string& vcf_file_name);

};


template <class ConsumerMT>
void VCFReaderMT<ConsumerMT>::parseFieldNames(const std::string& vcf_file_name) {

  std::ifstream vcf_file;

  // Open input file.

  vcf_file.open(vcf_file_name);

  if (not vcf_file.good()) {

    ExecEnv::log().critical("I/O error; could not open VCF file: {}", vcf_file_name);

  }

  try {

    long counter = 0;
    bool found_header = false;


    while (true) {

      std::string record_str;

      if (std::getline(vcf_file, record_str).eof()) break;

      std::string line_prefix = record_str.substr(0, FIELD_NAME_FRAGMENT_LENGTH_);

      if (line_prefix == FIELD_NAME_FRAGMENT_) {

        found_header = true;
        size_t field_count = 0;
        bt::char_separator<char> item_key_sep("\t");
        bt::tokenizer<bt::char_separator<char>> tokenize_item(record_str, item_key_sep);
        for(auto iter_item = tokenize_item.begin(); iter_item != tokenize_item.end(); ++iter_item) {

          if (field_count >= SKIP_FIELD_NAMES_) {

            field_names_.push_back(*iter_item);

          }

          ++field_count;

        }
        break;

      }

      ++counter;

    }

    vcf_file.close();

    if (not found_header) {

      ExecEnv::log().error("VCF Field Names Not Found");

    }

  }
  catch (std::exception const &e) {

    ExecEnv::log().critical("VCF file: {}, unexpected I/O exception: {}", vcf_file_name, e.what());

  }

}



template <class ConsumerMT>
const ContigId_t VCFReaderMT<ConsumerMT>::getContig(int32_t contig_idx) const {

  AutoMutex auto_mutex(mutex_); // lock on construction, unlock on destruction

  std::string contig_id = seqan::toCString(contigNames(context(*vcfIn_ptr_))[contig_idx]);

  return static_cast<const ContigId_t>(contig_id);

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
