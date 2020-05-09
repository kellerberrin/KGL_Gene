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

#include "kgl_variant_file_vcf_impl.h"


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////
// Dequeues VCF Records and passes them to the final parser logic which generates variant objects.

class VCFReaderMT {

public:

  explicit VCFReaderMT(const std::string& vcf_file_name) : vcf_io_(vcf_file_name) {}
  virtual ~VCFReaderMT() = default;

  // Perform multi-threaded parsing of queued VCF records.
  void readVCFFile();

  // Process each VCF record.
  virtual void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) = 0;

  // Process VCF header information.
  virtual void processVCFHeader(const VcfHeaderInfo& header_info) = 0;

  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return vcf_io_.getGenomeNames(); }

private:

  RecordVCFIO vcf_io_;                                   // VCF record queue.
  size_t consumer_thread_count_{8};                      // Consumer threads

  void readHeader();
  // Call the template VCF consumer class
  void VCFConsumer();

};



}   // end namespace


#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
