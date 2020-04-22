//
// Created by kellerberrin on 21/4/20.
//

#ifndef KGL_VARIANT_FILE_VCF_IMPL_H
#define KGL_VARIANT_FILE_VCF_IMPL_H


#include "kgl_genome_types.h"
#include "kgl_variant_file_impl.h"
#include "kgl_variant_file_vcf_record.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <memory>
#include <string>
#include <vector>

namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses the VCF header for Contig and Genome information.

class VCFParseHeader {

public:

  VCFParseHeader() = default;

  ~VCFParseHeader() = default;

  bool parseHeader(const std::string& file_name, std::unique_ptr<BaseStreamIO>&& vcf_stream);

  [[nodiscard]] const std::vector<std::string>& getGenomes() const { return vcf_genomes_; }

  [[nodiscard]] const VcfHeaderInfo& getHeaderInfo() const { return vcf_header_info_; }

private:

  std::vector<std::string> vcf_genomes_;                // Field (genome) names for each VCF record
  VcfHeaderInfo vcf_header_info_;
  // Parser constants.
  static constexpr const char* KEY_SEPARATOR_{"="};
  static constexpr const char* KEY_PREFIX_{"##"};
  static constexpr const char* FIELD_NAME_FRAGMENT_{"#CHROM"};
  static constexpr const size_t FIELD_NAME_FRAGMENT_LENGTH_{6};
  static constexpr const size_t SKIP_FIELD_NAMES_{9};  // Skip the fixed fields to the Genome names.


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF Record parser object (multi-threaded) - uses multiple threads to parse queued file lines from the IO reader.
// The line fields are tab delimited and are parsed into a VcfRecord which is then enqueued for further processing.

class RecordVCFIO : private FileVCFIO {

public:

  explicit RecordVCFIO(const std::string& vcf_file_name) : FileVCFIO(vcf_file_name),
                                                           vcf_record_queue_(HIGH_TIDE_, LOW_TIDE_) {}
  ~RecordVCFIO() override;

  void commenceVCFIO(size_t reader_threads) ; // Begin reading IO records, spawns threads.

  // export parser VCF records further up the parser chain.
  [[nodiscard]] std::unique_ptr<const VcfRecord> readVCFRecord() { return vcf_record_queue_.waitAndPop(); }

  // Read VCF header info.
  [[nodiscard]] const VcfHeaderInfo& VCFReadHeader() const;

  // Stored VCF header info.
  [[nodiscard]] const std::vector<std::string>& getGenomeNames() const { return parseheader_.getGenomes(); }

  [[nodiscard]] const std::string& getFileName() const { return fileName(); }

private:

  BoundedMtQueue<std::unique_ptr<const VcfRecord>> vcf_record_queue_; // Parsed VCF record queue
  std::vector<std::thread> vcf_record_thread_vec_;
  VCFParseHeader parseheader_;    // Get genome and contig information.
  size_t reader_threads_;

  static constexpr const long HIGH_TIDE_{10000};          // Maximum BoundedMtQueue size
  static constexpr const long LOW_TIDE_{1000};            // Low water mark to begin queueing VCF records
  static constexpr const long PARSER_THREADS_{4};         // Threads parsing into vcf_records.
  static constexpr const char HEADER_CHAR_{'#'};          // If first char start with '#' then a header record.
  static constexpr const size_t MINIMUM_VCF_FIELDS_{8};   // At least 8 fields, any others are format and genotype fields (header specified).
  static constexpr const char* VCF_FIELD_DELIMITER_{"\t"};   // VCF Field separator.
  const std::string FIELD_NOT_PRESENT_{"."}; // no field value

  void enqueueVCFRecord(); // enqueue vcf_records.
  bool parseVCFRecord(const std::unique_ptr<std::string>& line_record_ptr, const std::unique_ptr<VcfRecord>& vcf_record_ptr);
  bool moveToVcfRecord(std::vector<std::string>& fields, VcfRecord& vcf_record);

};



} // end namespace



#endif //KGL_KGL_VARIANT_FILE_VCF_IMPL_H
