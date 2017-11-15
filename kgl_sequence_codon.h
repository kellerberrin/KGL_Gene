//
// Created by kellerberrin on 14/11/17.
//

#ifndef KGL_SEQUENCE_CODON_H
#define KGL_SEQUENCE_CODON_H


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The codon class transfers indexed codons from a DNA5SequenceCoding + offset to the amino translation table
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <array>
#include "kgl_alphabet_base.h"
#include "kgl_sequence_base.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class Codon : public AlphabetSequence {

public:


  Codon(std::shared_ptr<DNA5SequenceCoding> sequence_ptr, ContigOffset_t codon_index);
  ~Codon() override = default;

  static ContigOffset_t codonLength(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) {

    return static_cast<ContigOffset_t>(sequence_ptr->length() / CODON_SIZE);

  }

  Nucleotide_DNA5_t operator[](size_t index) const { return bases_[index]; };

  std::string getSequenceAsString() const override {

    std::string codon_string;
    codon_string = bases_[0];
    codon_string += bases_[1];
    codon_string += bases_[2];
    return codon_string;

  }

  bool containsBaseN() const {

    return NucleotideColumn_DNA5::isNucleotideN(bases_[0])
           or NucleotideColumn_DNA5::isNucleotideN(bases_[1])
           or NucleotideColumn_DNA5::isNucleotideN(bases_[2]);

  }

  static constexpr ContigSize_t CODON_SIZE = 3;

private:

  std::array<Nucleotide_DNA5_t, CODON_SIZE> bases_;

};



} // namespace genome
} // namespace kellerberrin

#endif //KGL_SEQUENCE_CODON_H
