//
// Created by kellerberrin on 21/4/20.
//

#ifndef KGL_VARIANT_FILE_VCF_IMPL_H
#define KGL_VARIANT_FILE_VCF_IMPL_H


#include "kgl_genome_types.h"
#include "kgl_data_file_impl.h"
#include "kgl_variant_file_vcf_record.h"
#include "kgl_variant_factory_vcf_parse_header.h"

#include "kel_bound_queue.h"

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

using QueuedVCFRecord = std::optional<std::pair<size_t, std::unique_ptr<const VcfRecord>>>;

class RecordVCFIO {

public:

  RecordVCFIO() = default;
  ~RecordVCFIO() noexcept;

  // Begin reading IO records, spawns threads.
  bool commenceVCFIO(const std::string& vcf_file_name);

  // export parser VCF records further up the parser chain.
  [[nodiscard]] QueuedVCFRecord readVCFRecord() { return vcf_record_queue_.waitAndPop(); }

  // Push an eof marker onto the queue
  void enqueueEOF() { vcf_record_queue_.push(QUEUED_EOF_MARKER); }

  [[nodiscard]] [[maybe_unused]] const std::string& getFileName() const { return file_data_.fileName(); }


private:

  // The downstream queue of line records.
  FileDataIO file_data_;

  // VCF queue parameters.
  static constexpr const size_t VCF_HIGH_TIDE_{2000};     // Maximum BoundedMtQueue size
  static constexpr const size_t VCF_LOW_TIDE_{1000};       // Low water mark to begin queueing VCF records
  static constexpr const char* VCF_NAME_{"VCF Record Queue"};    // The queue name
  static constexpr const size_t VCF_SAMPLE_RATE_{500}; // Queue sample rate.
  // Parsed VCF record queue
  BoundedMtQueue<QueuedVCFRecord> vcf_record_queue_{VCF_HIGH_TIDE_, VCF_LOW_TIDE_, VCF_NAME_, VCF_SAMPLE_RATE_};

  // VCF queue worker threads
  static constexpr const long PARSER_THREADS_{15};         // Threads parsing into vcf_records.
  // The detached main thread.
  ThreadPool detached_launch_{1};
  // Synchronize shutdown
  std::future<void> launch_token_;

  // VCF record constants.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record (skip).
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char VCF_FIELD_DELIMITER_CHAR_{'\t'};   // VCF Field separator (char).
  const std::string FIELD_NOT_PRESENT_{"."}; // no field value

  void launchThreads();
  void enqueueVCFRecord(); // enqueue vcf_records.
  std::optional<std::unique_ptr<VcfRecord>> moveToVcfRecord(std::string&& line_record);

};



} // end namespace



#endif //KGL_KGL_VARIANT_FILE_VCF_IMPL_H
