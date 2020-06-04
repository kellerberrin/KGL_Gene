//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_genome_types.h"
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
                   InfoDataEvidence info_data_block,
                   std::optional<std::shared_ptr<FormatData>> format_data,
                   uint32_t alternate_variant_index = 0,
                   uint32_t alternate_variant_count = 1)
                   : vcf_record_count_(vcf_record_count),
                     info_data_block_(std::move(info_data_block)),
                     format_data_(std::move(format_data)),
                     alternate_variant_index_(alternate_variant_index),
                     alternate_variant_count_(alternate_variant_count) {}

  VariantEvidence(const VariantEvidence&) = default;
  ~VariantEvidence() {

     if (format_data_) {

       format_data_ = std::nullopt;

     }

     if (info_data_block_) {

       info_data_block_ = std::nullopt;

     }

  }

  [[nodiscard]] std::string output(char delimiter, VariantOutputIndex output_index) const;

  [[nodiscard]] size_t vcfRecordCount() const { return vcf_record_count_; }

  [[nodiscard]] InfoDataEvidence infoData() const { return info_data_block_; }

  [[nodiscard]] std::optional<std::shared_ptr<const FormatData>> formatData() const { return format_data_; }

  [[nodiscard]] uint32_t altVariantIndex() const { return alternate_variant_index_; }

  [[nodiscard]] uint32_t altVariantCount() const { return alternate_variant_count_; }

private:

  size_t vcf_record_count_; // The VCF line count, the original file line record that generated this variant.
  InfoDataEvidence info_data_block_;   // INFO data items, may be missing.
  std::optional<std::shared_ptr<FormatData>> format_data_;  // Format data items, may be missing.
  // Zero based index. Which of the alternate variants (from left to right in the comma delimited alt field) is this variant.
  // These variables can be used to access Info field vectors that are designated Type='A' for alternate allele.
  uint32_t alternate_variant_index_; // The default index 0 / count 1 implies 1 alternate variant (the usual case).
  uint32_t alternate_variant_count_; // How many comma delimited alternate variants were specified in the VCF record.

};




}   // end namespace


#endif //KGL_VARIANT_EVIDENCE_H
