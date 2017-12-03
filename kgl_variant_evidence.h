//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_variant.h"



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

  virtual bool isReadCount() const { return false; }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The variant read count evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ReadCountEvidence : public VariantEvidence {

public:

  ReadCountEvidence(NucleotideReadCount_t read_count,
                    NucleotideReadCount_t mutant_count,
                    NucleotideReadCount_t const count_array[],
                    ContigSize_t count_array_size) : read_count_(read_count),
                                                     mutant_count_(mutant_count) {

    for(size_t idx = 0; idx < count_array_size; ++idx) {

      count_array_.push_back(count_array[idx]);

    }

  }
  ReadCountEvidence(const ReadCountEvidence& variant) = default;
  ~ReadCountEvidence() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return mutant_count_; }
  const std::vector<NucleotideReadCount_t>& countArray() const { return count_array_; }

  double proportion() const { return static_cast<double>(mutant_count_) / static_cast<double>(read_count_); }

  std::string output(char delimiter, VariantOutputIndex output_index) const override;

  bool isReadCount() const override { return true; }

private:

  NucleotideReadCount_t read_count_;
  NucleotideReadCount_t mutant_count_;
  std::vector<NucleotideReadCount_t> count_array_;

};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_EVIDENCE_H
