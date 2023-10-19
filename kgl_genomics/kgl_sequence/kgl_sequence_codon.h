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


namespace kellerberrin::genome {   //  organization level namespace


class Codon  {

public:


  Codon(const DNA5SequenceCoding& coding_sequence, ContigOffset_t codon_index);
  ~Codon() = default;

  [[nodiscard]] static ContigSize_t codonLength(const DNA5SequenceCoding& coding_sequence) {

    return static_cast<ContigSize_t>(coding_sequence.length() / CODON_SIZE);

  }

  [[nodiscard]] CodingDNA5::Alphabet operator[](size_t index) const { return bases_[index]; };

  [[nodiscard]] std::string getSequenceAsString() const {

    std::string codon_string;
    codon_string = CodingDNA5::convertToChar(bases_[0]);
    codon_string += CodingDNA5::convertToChar(bases_[1]);
    codon_string += CodingDNA5::convertToChar(bases_[2]);
    return codon_string;

  }

  [[nodiscard]] bool containsBaseN() const {

    return bases_[0] == CodingDNA5::Alphabet::N
           or bases_[1] == CodingDNA5::Alphabet::N
           or bases_[2]  == CodingDNA5::Alphabet::N;

  }

  static constexpr ContigSize_t CODON_SIZE = 3;

private:

  std::array<CodingDNA5::Alphabet, CODON_SIZE> bases_;

};



} // namespace genome


#endif //KGL_SEQUENCE_CODON_H
