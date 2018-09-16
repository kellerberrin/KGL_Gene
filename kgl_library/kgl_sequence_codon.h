//
// Created by kellerberrin on 14/11/17.
//

#ifndef KGL_SEQUENCE_CODON_H
#define KGL_SEQUENCE_CODON_H


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The codon class transfers indexed codons from a DNA5SequenceCoding + offset to the amino translation table
// Only accepts the stranded sequence DNA5SequenceCoding.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <array>
#include "kgl_sequence_base.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class Codon  {

public:


  Codon(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr, ContigOffset_t codon_index);
  ~Codon() = default;

  static ContigSize_t codonLength(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) {

    return static_cast<ContigSize_t>(sequence_ptr->length() / CODON_SIZE);

  }

  CodingDNA5::Alphabet operator[](size_t index) const { return bases_[index]; };
  void modifyBase(size_t index, CodingDNA5::Alphabet base);

  std::string getSequenceAsString() const {

    std::string codon_string;
    codon_string = CodingDNA5::convertToChar(bases_[0]);
    codon_string += CodingDNA5::convertToChar(bases_[1]);
    codon_string += CodingDNA5::convertToChar(bases_[2]);
    return codon_string;

  }

  bool containsBaseN() const {

    return bases_[0] == CodingDNA5::Alphabet::N
           or bases_[1] == CodingDNA5::Alphabet::N
           or bases_[2]  == CodingDNA5::Alphabet::N;

  }

  static constexpr ContigSize_t CODON_SIZE = 3;

private:

  std::array<CodingDNA5::Alphabet, CODON_SIZE> bases_;

};



} // namespace genome
} // namespace kellerberrin


#endif //KGL_SEQUENCE_CODON_H
