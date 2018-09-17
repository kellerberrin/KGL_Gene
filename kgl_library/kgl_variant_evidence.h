//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include <sstream>
#include "kgl_genome_types.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This class holds the evidence that resulted in the creation of a variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantEvidence { // Top level object.

public:

  explicit VariantEvidence() = default;
  virtual ~VariantEvidence() = default;

  virtual std::string output(char delimiter, VariantOutputIndex output_index) const = 0;

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Basic count evidence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class CountEvidence : public VariantEvidence {

public:

  explicit CountEvidence(size_t ref_count,
                         size_t alt_count,
                         size_t DP_count,
                         Phred_t GQ_value,
                         Phred_t quality,
                         size_t vcf_record_count) : ref_count_(ref_count),
                                                    alt_count_(alt_count),
                                                    DP_count_(DP_count),
                                                    GQ_value_(GQ_value),
                                                    quality_(quality),
                                                    vcf_record_count_(vcf_record_count) {}
  virtual ~CountEvidence() = default;

  size_t refCount() const { return ref_count_; }
  size_t altCount() const { return alt_count_; }
  size_t DPCount() const { return DP_count_; }
  Phred_t GQValue() const { return GQ_value_; }
  Phred_t Quality() const { return quality_; }
  size_t vcfRecordCount() const { return vcf_record_count_; }

  std::string output(char delimiter, VariantOutputIndex) const override;

private:

  size_t ref_count_;
  size_t alt_count_;
  size_t DP_count_;
  Phred_t GQ_value_;
  Phred_t quality_;
  size_t vcf_record_count_;

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_EVIDENCE_H
