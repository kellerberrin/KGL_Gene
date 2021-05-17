//
// Created by kellerberrin on 28/02/18.
//

#ifndef KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H


#include "kel_utility.h"
#include "kgl_resource_db.h"
#include "kgl_variant_file_vcf_record.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF record parser.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParseVCFRecord {

public:

  ParseVCFRecord(const VcfRecord& vcf_record, const std::shared_ptr<const GenomeReference>& genome_db_ptr);
  ~ParseVCFRecord() = default;

  [[nodiscard]] const std::vector<std::string>& formatFields() const { return format_fields_; }
  [[nodiscard]] Phred_t quality() const { return quality_; }
  [[nodiscard]] const std::string& reference() const { return reference_; }
  [[nodiscard]] const std::vector<std::string>& alleles() const { return alleles_; }
  [[nodiscard]] ContigOffset_t offset() const { return allele_offset_; }
  [[nodiscard]] const std::shared_ptr<const ContigReference>& contigPtr() const { return contig_ptr_; }
  [[nodiscard]] bool passedFilter() const { return passed_filter_; }
  [[nodiscard]] bool parseResult() const { return parse_result_; }

  [[nodiscard]] std::optional<size_t> formatIndex(const std::string& format_code) const;

  constexpr static const char* FORMAT_GT{"GT"};
  constexpr static const char* FORMAT_GQ{"GQ"};
  constexpr static const char* FORMAT_PL{"PL"};    // Gatk uses PL (integers)
  constexpr static const char* FORMAT_GL{"GL"};    // FreeBayes uses GL (floats)
  constexpr static const char* FORMAT_AD{"AD"};
  constexpr static const char* FORMAT_DP{"DP"};

private:

  std::vector<std::string> format_fields_;  // GT:AD:DP:GQ:PGT:PID:PL
  std::string reference_;
  std::vector<std::string> alleles_;
  ContigOffset_t allele_offset_{0};
  std::shared_ptr<const ContigReference> contig_ptr_;
  Phred_t quality_{0.0};
  bool passed_filter_{false};
  bool parse_result_{true};

  constexpr static const char FORMAT_SEPARATOR_{':'};
  constexpr static const char ALLELE_SEPARATOR_{','};
  constexpr static const char* PASSED_FILTERS_{"PASS"};

};




}   // end namespace



#endif //KGL_KGL_VARIANT_FACTORY_RECORD_VCF_IMPL_H
