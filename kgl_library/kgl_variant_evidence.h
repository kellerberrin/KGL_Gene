//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include "kgl_genome_types.h"

#include <sstream>
#include <string>
#include <memory>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class holds the evidence that resulted in the creation of a variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantEvidence { // Top level object.

public:

  explicit VariantEvidence(std::shared_ptr<std::string> info) : info_(std::move(info)) {}
  virtual ~VariantEvidence() = default;

  [[nodiscard]] virtual std::string output(char delimiter, VariantOutputIndex) const { return delimiter  + *info_ + delimiter; }

  [[nodiscard]] const std::string& info_field() const { return *info_; }

private:

  // The VCF file info field.
  std::shared_ptr<std::string> info_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic count evidence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CountEvidence : public VariantEvidence {

public:

  explicit CountEvidence(std::shared_ptr<std::string> info,
                         size_t ref_count,
                         size_t alt_count,
                         size_t DP_count,
                         Phred_t GQ_value,
                         Phred_t quality,
                         size_t vcf_record_count) : VariantEvidence(std::move(info)),
                                                    ref_count_(ref_count),
                                                    alt_count_(alt_count),
                                                    DP_count_(DP_count),
                                                    GQ_value_(GQ_value),
                                                    quality_(quality),
                                                    vcf_record_count_(vcf_record_count) {}
  ~CountEvidence() override = default;

  [[nodiscard]] size_t refCount() const { return ref_count_; }
  [[nodiscard]] size_t altCount() const { return alt_count_; }
  [[nodiscard]] size_t DPCount() const { return DP_count_; }
  [[nodiscard]] Phred_t GQValue() const { return GQ_value_; }
  [[nodiscard]] Phred_t Quality() const { return quality_; }
  [[nodiscard]] size_t vcfRecordCount() const { return vcf_record_count_; }

  [[nodiscard]] std::string output(char delimiter, VariantOutputIndex) const override;

private:

  size_t ref_count_;
  size_t alt_count_;
  size_t DP_count_;
  Phred_t GQ_value_;
  Phred_t quality_;
  size_t vcf_record_count_;

};


}   // end namespace


#endif //KGL_VARIANT_EVIDENCE_H
