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
#include "kel_queue_mt_safe.h"
#include "kel_workflow_threads.h"

#include "kgl_variant_vcf_impl.h"


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
  virtual void ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) = 0;

  // Process VCF header information.
  virtual void processVCFHeader(const VCFHeaderInfo& header_info) = 0;

  // Stored VCF header info.
  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  constexpr static const size_t DEFAULT_PARSER_THREADS{50};

private:

  // VCF record queue.
  ParseVCF vcf_parser_;

  // Threads to process the VCF record queue.
  WorkflowThreads parser_threads_;

  // Get genome and contig_ref_ptr information.
  VCFParseHeader parseheader_;

  void readHeader(const std::string& file_name);
  // Call the template VCF consumer class
  void VCFConsumer();

};



}   // end namespace


#endif //KGL_VARIANT_FACTORY_READVCF_IMPL_H
