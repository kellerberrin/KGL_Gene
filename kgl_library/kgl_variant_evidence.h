//
// Created by kellerberrin on 4/12/17.
//

#ifndef KGL_VARIANT_EVIDENCE_H
#define KGL_VARIANT_EVIDENCE_H


#include <sstream>
#include "kgl_genome_types.h"
#include "kgl_alphabet_readcount.h"



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

  Phred_t Quality() const { return 100.0; }

  bool isReadCount() const { return true; }

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

  Phred_t Quality() const { return quality_; }
  bool isVCF() const override { return true; }

private:

  std::string info_;
  Phred_t quality_;

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
