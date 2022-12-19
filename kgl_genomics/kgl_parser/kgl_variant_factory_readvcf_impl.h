//
// Created by kellerberrin on 26/02/18.
//

#ifndef KGL_VARIANT_FACTORY_READVCF_IMPL_H
#define KGL_VARIANT_FACTORY_READVCF_IMPL_H


#include <string>
#include <thread>
#include <fstream>
#include <functional>

#include "kel_exec_env.h"
#include "../../kel_thread/kel_mt_queue.h"
#include "../../kel_thread/kel_thread_pool.h"

#include "kgl_variant_file_vcf_impl.h"


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////
// Dequeues VCF Records and passes them to the final parser logic which generates variant objects.

class VCFReaderMT {

public:

  explicit VCFReaderMT(size_t thread_count = DEFAULT_PARSER_THREADS) : parser_threads_(thread_count) {}
  virtual ~VCFReaderMT() = default;

  // Perform multi-threaded parsing of queued VCF records.
  void readVCFFile(const std::string& vcf_file_name);

  // Process each VCF record.
  virtual void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) = 0;

  // Process VCF header information.
  virtual void processVCFHeader(const VcfHeaderInfo& header_info) = 0;

  // Stored VCF header info.
  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  constexpr static const size_t DEFAULT_PARSER_THREADS{50};

private:

  // VCF record queue.
  RecordVCFIO vcf_io_;

  // Threads to process the VCF record queue.
  WorkflowThreads parser_threads_;

  // Get genome and contig information.
  VCFParseHeader parseheader_;

  void readHeader(const std::string& file_name);
  // Call the template VCF consumer class
  void VCFConsumer();

};



}   // end namespace


#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
