//
// Created by kellerberrin on 14/11/17.
//

#ifndef KGL_SEQUENCE_CODON_H
#define KGL_SEQUENCE_CODON_H


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The codon class transfers indexed codons from a DNA5SequenceCoding + offset to the amino translation table
// Only accepts the stranded sequence DNA5SequenceCoding.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <array>
#include <format>
#include "kgl_sequence_base.h"


namespace kellerberrin::genome {   //  organization level namespace


class Codon  {

public:


  Codon(const DNA5SequenceCoding& coding_sequence, ContigOffset_t codon_index);
  ~Codon() = default;

  /// ReturnType the number of whole codons in the coding sequence.
  [[nodiscard]] static ContigSize_t codonLength(const DNA5SequenceCoding& coding_sequence) {

    return static_cast<ContigSize_t>(coding_sequence.length() / CODON_SIZE);

  }

  /// Random access to the codon bases. Bounds checked: an out of range index throws std::out_of_range.
  [[nodiscard]] CodingDNA5::Alphabet operator[](size_t index) const { return bases_.at(index); }

  /// ReturnType the codon as a 3 character base string.
  [[nodiscard]] std::string getSequenceAsString() const {

    return std::format("{}{}{}",
                       CodingDNA5::convertToChar(bases_[0]),
                       CodingDNA5::convertToChar(bases_[1]),
                       CodingDNA5::convertToChar(bases_[2]));

  }

  /// True if the codon contains the unknown base 'N'.
  [[nodiscard]] bool containsBaseN() const {

    return bases_[0] == CodingDNA5::Alphabet::N
           or bases_[1] == CodingDNA5::Alphabet::N
           or bases_[2] == CodingDNA5::Alphabet::N;

  }

  /// The codon size in bases.
  static constexpr ContigSize_t CODON_SIZE = 3;
  /// ReturnType the remainder bases of a sequence size that are not part of a whole codon.
  [[nodiscard]] static size_t codonRemainder(ContigSize_t size) { return size % CODON_SIZE; }

private:

  std::array<CodingDNA5::Alphabet, CODON_SIZE> bases_;

};



} // namespace genome


#endif //KGL_SEQUENCE_CODON_H