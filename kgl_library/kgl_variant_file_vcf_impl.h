//
// Created by kellerberrin on 21/4/20.
//

#ifndef KGL_VARIANT_FILE_VCF_IMPL_H
#define KGL_VARIANT_FILE_VCF_IMPL_H


#include "kgl_genome_types.h"
#include "kgl_variant_file_impl.h"
#include "kgl_variant_file_vcf_record.h"
#include "kgl_variant_factory_vcf_parse_header.h"

#include <memory>
#include <string>
#include <vector>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF Record parser object (multi-threaded) - uses multiple threads to parse queued file lines from the IO reader.
// The line fields are tab delimited and are parsed into a VcfRecord which is then enqueued for further processing.

using QueuedVCFRecord = std::optional<std::pair<size_t, std::unique_ptr<const VcfRecord>>>;

class RecordVCFIO : private FileVCFIO {

public:

  explicit RecordVCFIO(const std::string& vcf_file_name) : FileVCFIO(vcf_file_name),
                                                           vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~RecordVCFIO() override;

  void commenceVCFIO(size_t reader_threads) ; // Begin reading IO records, spawns threads.

  // export parser VCF records further up the parser chain.
  [[nodiscard]] QueuedVCFRecord readVCFRecord() { return vcf_record_queue_.waitAndPop(); }

  [[nodiscard]] const std::string& getFileName() const { return fileName(); }

private:

  BoundedMtQueue<QueuedVCFRecord> vcf_record_queue_; // Parsed VCF record queue
  std::vector<std::thread> vcf_record_thread_vec_;
  size_t reader_threads_;

  static constexpr const long HIGH_TIDE_{100000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{10000};            // Low water mark to begin queueing VCF records
  static constexpr const long PARSER_THREADS_{4};         // Threads parsing into vcf_records.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record (skip).
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char* VCF_FIELD_DELIMITER_{"\t"};   // VCF Field separator.
  static constexpr const char VCF_FIELD_DELIMITER_CHAR_{'\t'};   // VCF Field separator (char).
  const std::string FIELD_NOT_PRESENT_{"."}; // no field value

  void enqueueVCFRecord(); // enqueue vcf_records.
  bool parseVCFRecord(std::unique_ptr<const std::string> line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr);
  bool moveToVcfRecord(std::unique_ptr<const std::string> line_record_ptr, std::vector<std::string>& fields, VcfRecord& vcf_record);

};



} // end namespace



#endif //KGL_KGL_VARIANT_FILE_VCF_IMPL_H
