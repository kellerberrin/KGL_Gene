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
  virtual Phred_t calculateQuality() const = 0;

  virtual bool isReadCount() const { return false; }
  virtual bool isVCF() const { return false; }

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The variant read count evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class Alphabet>
class ReadCountEvidence : public VariantEvidence {

public:

  ReadCountEvidence(size_t mutant_offset,
                    NucleotideReadCount_t const count_array[]) : mutant_offset_(mutant_offset) {

    // check mutant offset.
    if (mutant_offset >= Alphabet::NUCLEOTIDE_COLUMNS) {

      ExecEnv::log().error("ReadCountEvidence; Bad mutant index: {} for count_array size: {}", mutant_offset, Alphabet::NUCLEOTIDE_COLUMNS);
      mutant_offset = 0;

    }
    read_count_ = 0;
    for(size_t idx = 0; idx < Alphabet::NUCLEOTIDE_COLUMNS; ++idx) {

      read_count_ += count_array[idx];
      count_array_.push_back(count_array[idx]);

    }

  }
  ReadCountEvidence(const ReadCountEvidence& variant) = default;
  ~ReadCountEvidence() override = default;

  NucleotideReadCount_t readCount() const { return read_count_; }
  NucleotideReadCount_t mutantCount() const { return count_array_[mutant_offset_]; }
  const std::vector<NucleotideReadCount_t>& countArray() const { return count_array_; }

  double proportion() const { return static_cast<double>(mutantCount()) / static_cast<double>(readCount()); }

  std::string output(char delimiter, VariantOutputIndex output_index) const override;

  virtual Phred_t calculateQuality() const override { return 100.0; }

  bool isReadCount() const override { return true; }

private:

  NucleotideReadCount_t read_count_;
  size_t mutant_offset_;
  std::vector<NucleotideReadCount_t> count_array_;

};


template <class Alphabet>
std::string ReadCountEvidence<Alphabet>::output(char delimiter, VariantOutputIndex) const
{

  std::stringstream ss;

  ss << mutantCount() << "/" << readCount() << delimiter;
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << ReadCountColumns::convertToChar(ReadCountColumns::offsetToNucleotide(idx)) << ":" << countArray()[idx] << delimiter;
  }

  return ss.str();

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The VCF evidence object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VCFEvidence : public VariantEvidence {

public:

  explicit VCFEvidence(const std::string& info, Phred_t quality) : info_(info), quality_(quality) {}
  virtual ~VCFEvidence() = default;

  std::string output(char, VariantOutputIndex) const override { return info_; }
  Phred_t calculateQuality() const override { return quality_; }
  bool isVCF() const override { return true; }

private:

  std::string info_;
  Phred_t quality_;
  Haplotypes haplotypes_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_EVIDENCE_H
