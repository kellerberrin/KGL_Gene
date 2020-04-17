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


namespace kellerberrin::genome {   //  organization::project level namespace



class VCFReaderMT {

public:

  explicit VCFReaderMT(const std::string& vcf_file_name) : producer_consumer_queue_(HIGH_TIDE_, LOW_TIDE_),
                                                           vcf_io_ptr_(std::make_unique<SeqanVCFIO>(vcf_file_name)) {}

  virtual ~VCFReaderMT() = default;

  // Perform multi-threaded parsing of the VCF file.
  void readVCFFile();

  // Process each VCF record.
  virtual void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) = 0;

  // Process VCF header information.
  virtual void processVCFHeader(const VcfHeaderInfo& header_info) = 0;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return vcf_io_ptr_->getGenomeNames(); }

private:

  // The bounded queue data type.
  using ReaderQueueRecord = std::pair<size_t, std::unique_ptr<VcfRecord>>;

  BoundedMtQueue<ReaderQueueRecord> producer_consumer_queue_; // The Producer/Consumer record queue
  std::unique_ptr<BaseVCFIO> vcf_io_ptr_;

  size_t consumer_thread_count_{2};                      // Consumer threads (defaults to local CPU cores available or max)
  static constexpr const int MAX_CONSUMER_THREADS_{4};     // Max consumer threads. Spawning more threads does not increase performance
  static constexpr const int MIN_CONSUMER_THREADS_{1};     // Need at least 1 consumer thread

  static constexpr const size_t report_increment_{10000};

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records


  void readHeader();
  // Read the VCF file and queue the record in a BoundedMtQueue.
  void VCFProducer();
  // Call the template VCF consumer class
  void VCFConsumer();


};



}   // end namespace


#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
