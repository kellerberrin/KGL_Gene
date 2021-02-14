//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_genome_types.h"
#include "kgl_variant_db_type.h"
#include "kgl_variant_factory_vcf_evidence.h"

#include <sstream>
#include <string>
#include <memory>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class holds the evidence that resulted in the creation of a variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Per Variant format data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class FormatData  {

public:

   FormatData( size_t ref_count,
               size_t alt_count,
               size_t DP_count,
               Phred_t GQ_value,
               Phred_t quality) : ref_count_(ref_count),
                                  alt_count_(alt_count),
                                  DP_count_(DP_count),
                                  GQ_value_(GQ_value),
                                  quality_(quality) {}
  ~FormatData() = default;

  [[nodiscard]] size_t refCount() const { return ref_count_; }
  [[nodiscard]] size_t altCount() const { return alt_count_; }
  [[nodiscard]] size_t DPCount() const { return DP_count_; }
  [[nodiscard]] Phred_t GQValue() const { return GQ_value_; }
  [[nodiscard]] Phred_t Quality() const { return quality_; }

  [[nodiscard]] std::string output(char delimiter) const;

private:

  size_t ref_count_;
  size_t alt_count_;
  size_t DP_count_;
  Phred_t GQ_value_;
  Phred_t quality_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VariantEvidence { // Top level object.

public:

  VariantEvidence( size_t vcf_record_count,
                   DataSourceEnum data_source,
                   bool pass_filter,
                   std::shared_ptr<const DataMemoryBlock> info_data_block,
                   std::shared_ptr<const FormatData> format_data,
                   uint32_t alternate_variant_index = 0,
                   uint32_t alternate_variant_count = 1)
                   : vcf_record_count_(vcf_record_count),
                     data_source_(data_source),
                     pass_filter_(pass_filter),
                     info_data_block_(std::move(info_data_block)),
                     format_data_(std::move(format_data)),
                     alternate_variant_index_(alternate_variant_index),
                     alternate_variant_count_(alternate_variant_count) {}

  VariantEvidence() = default;
  VariantEvidence(const VariantEvidence& copy) = default;
  ~VariantEvidence() = default;

  [[nodiscard]] std::string output(char delimiter, VariantOutputIndex output_index) const;

  [[nodiscard]] size_t vcfRecordCount() const { return vcf_record_count_; }

  [[nodiscard]] bool passFilter() const { return pass_filter_; }

  [[nodiscard]] DataSourceEnum dataSource() const { return data_source_; }

  // Must be checked for null pointer before use.
  [[nodiscard]] std::shared_ptr<const DataMemoryBlock> infoData() const { return info_data_block_; }

  // Must be checked for null pointer before use.
  [[nodiscard]] std::shared_ptr<const FormatData> formatData() const { return format_data_; }

  [[nodiscard]] uint32_t altVariantIndex() const { return alternate_variant_index_; }

  [[nodiscard]] uint32_t altVariantCount() const { return alternate_variant_count_; }


private:

  size_t vcf_record_count_{0}; // The VCF line count, the original file line record that generated this variant.
  DataSourceEnum data_source_{DataSourceEnum::NotImplemented};
  bool pass_filter_{true};  // The VCF record has "PASS"ed all quality filters.
  std::shared_ptr<const DataMemoryBlock> info_data_block_;   // INFO data items, may be missing.
  std::shared_ptr<const FormatData> format_data_;  // Format data items, may be missing.
  // Zero based index. Which of the alternate variants (from left to right in the comma delimited alt field) is this variant.
  // These variables can be used to access Info field vectors that are designated Type='A' for alternate allele.
  uint32_t alternate_variant_index_{0}; // The default index 0 / count 1 implies 1 alternate variant (the usual case).
  uint32_t alternate_variant_count_{0}; // How many comma delimited alternate variants were specified in the VCF record.


};




}   // end namespace


#endif //KGL_VARIANT_EVIDENCE_H
