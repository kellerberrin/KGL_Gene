//
// Created by kellerberrin on 20/10/17.
//

#ifndef SAMFILE_KGL_SEQUENCE_H
#define SAMFILE_KGL_SEQUENCE_H

#include <string>
#include "kgl_nucleotide.h"
#include "kgl_base_sequence.h"
#include "kgl_amino.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
class AminoSequence {

public:

  using AminoType = typename T::AminoType;
  using ProteinString = std::basic_string<AminoType>;

  explicit AminoSequence(ProteinString sequence) : amino_sequence_(std::move(sequence)) {};
  AminoSequence() = delete;
  virtual ~AminoSequence() = default;

  AminoType operator[] (ContigOffset_t& offset) const { return amino_sequence_[offset]; }
  ContigSize_t length() const { return amino_sequence_.length(); }
  const AminoType* baseAddress(ContigOffset_t& offset) const { return &amino_sequence_[offset]; }


private:

  ProteinString amino_sequence_;

};

using StandardAminoSequence = AminoSequence<AminoAcidTypes>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Coding Sequence - An assembly of DNA/RNA sequences from CDS features, should contain start and stop codons.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class T, class AminoTable>
class CodingSequence {

public:


  explicit CodingSequence(std::shared_ptr<BaseSequence<T>> sequence_ptr) : sequence_ptr_(sequence_ptr) {};
  CodingSequence() = delete;
  ~CodingSequence() = default;

  static std::string translationTableName() { return AminoTable::TableName(); }
  inline AminoAcidTypes::Codon firstCodon() const { return getCodon(0); }
  bool checkStartCodon() const { return AminoTable::isStartCodon(firstCodon()); }
  inline AminoAcidTypes::Codon lastCodon() const { return getCodon(codonLength() - 1); }
  bool checkStopCodon() const { return AminoTable::isStopCodon(lastCodon()); }
  size_t checkNonsenseMutation() const {

    for (size_t index = 0; index < codonLength() - 1; ++index) {

      if (AminoTable::isStopCodon(getCodon(index))) return index;

    }

    return 0;

  }
  inline size_t codonLength() const { return sequence_ptr_->length() / 3; }
  inline AminoAcidTypes::Codon getCodon(size_t index) const {

    if (index >= codonLength()) {

      ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                           index, sequence_ptr_->length());
      index = 0;
    }

    AminoAcidTypes::Codon codon;
    index = index * 3;
    codon.bases = sequence_ptr_->baseAddress(index);
    return codon;

  }

private:

  std::shared_ptr<BaseSequence<T>> sequence_ptr_;

};


using DNA5Sequence = BaseSequence<NucleotideColumn_DNA5>;
using StandardCodingSequence = CodingSequence<NucleotideColumn_DNA5, StandardTranslationTable_1>;


}   // namespace genome
}   // namespace kellerberrin

#endif //KGL_SEQUENCE_H
