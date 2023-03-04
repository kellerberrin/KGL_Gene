//
// Created by kellerberrin on 3/03/23.
//

#ifndef KGL_VARIANT_VCF_IMPL_H
#define KGL_VARIANT_VCF_IMPL_H


#include "kgl_genome_types.h"
#include "kel_basic_io.h"
#include "kgl_variant_vcf_record.h"
#include "kgl_variant_factory_vcf_parse_header.h"

#include "kel_queue_tidal.h"
#include "kel_workflow_pipeline.h"

#include <memory>
#include <string>
#include <vector>

namespace kellerberrin::genome {   //  organization::project level namespace


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// VCF Record parser object (multi-threaded) - uses multiple threads to parse queued file lines from the IO reader.
// The line fields are tab delimited and are parsed into a VcfRecord which is then enqueued for further processing.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using VCFRecordPtr = std::unique_ptr<const VcfRecord>;

class ParseVCF {

  using VCFPipeline = WorkflowPipeline<IOLineRecord, VCFRecordPtr>;

public:

  ParseVCF() = default;
  ~ParseVCF();

  // Begin reading IO records, allocates threads to bgz decompression and parsing the VCF record.
  bool open( const std::string& vcf_file_name,
             size_t decompression_threads = DECOMPRESSION_THREADS_,
             size_t vcf_parse_threads = PARSER_THREADS_);

  // Export parsed VCF records further up the parser chain.
  // Note that nullptr indicates an EOF condition.
  [[nodiscard]] VCFRecordPtr readVCFRecord() { return vcf_pipeline_.waitAndPop(); }

  // Push an eof marker onto the queue
  void enqueueEOF() { vcf_pipeline_.push(IOLineRecord::createEOFMarker()); }

  [[nodiscard]] const std::string& getFileName() const { return file_name_; }

private:

  // Stream to read decompressed VCF records.
  std::unique_ptr<BaseStreamIO> stream_ptr_;
  std::string file_name_;
  // VCF queue worker threads
  static constexpr const long DECOMPRESSION_THREADS_{15};         // Threads decompressing bgz records.
  // VCF queue worker threads
  static constexpr const long PARSER_THREADS_{15};         // Threads parsing into vcf_records.
  // Pipeline to parse the VCF line record into fields.
  VCFPipeline vcf_pipeline_;
  // Read the decompressed line records and enqueue them in the VCF parser pipeline (1 thread).
  WorkflowThreads enqueue_thread_{1};

  // VCF record constants.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record (skip).
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char VCF_FIELD_DELIMITER_CHAR_{'\t'};   // VCF Field separator (char).
  const std::string FIELD_NOT_PRESENT_{"."}; // no field value
  // Mandatory Field Offsets.
  static constexpr const size_t CONTIG_FIELD_IDX_{0};
  static constexpr const size_t OFFSET_FIELD_IDX_{1};
  static constexpr const size_t IDENT_FIELD_IDX_{2};
  static constexpr const size_t REF_FIELD_IDX_{3};
  static constexpr const size_t ALT_FIELD_IDX_{4};
  static constexpr const size_t QUALITY_FIELD_IDX_{5};
  static constexpr const size_t FILTER_FIELD_IDX_{6};
  static constexpr const size_t INFO_FIELD_IDX_{7};
  static constexpr const size_t FORMAT_FIELD_IDX_{8};

  void enqueueLineRecord();
  std::unique_ptr<VcfRecord> moveToVcfRecord(IOLineRecord line_record);

};



} // end namespace


#endif //KGL_VARIANT_VCF_IMPL_H
