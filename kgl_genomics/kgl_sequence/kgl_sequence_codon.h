//
// kgl_sequence_codon.h — a plain value type over three Nucleotides.
//

#ifndef KGL_SEQUENCE_CODON_H
#define KGL_SEQUENCE_CODON_H


#include <array>
#include <optional>
#include <span>
#include <string>

#include "kgl_genome_types.h"
#include "kgl_alphabet.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class DNA5SequenceCoding;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Codon — transfers an indexed codon from a DNA5SequenceCoding to the amino translation table.
// Only accepts the stranded DNA5SequenceCoding sequence.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Codon {

public:

  Codon() = default;

  /// The codon size in bases.
  static constexpr ContigSize_t CODON_SIZE = 3;

  /// ReturnType the number of whole codons in the coding sequence.
  [[nodiscard]] static ContigSize_t codonLength(ContigSize_t coding_sequence_size) noexcept {
    return static_cast<ContigSize_t>(coding_sequence_size / CODON_SIZE);
  }
  /// ReturnType the number of whole codons in the coding sequence (retained sequence overload).
  [[nodiscard]] static ContigSize_t codonLength(const DNA5SequenceCoding& coding_sequence) noexcept;
  /// ReturnType the remainder bases of a sequence size that are not part of a whole codon.
  [[nodiscard]] static std::size_t codonRemainder(ContigSize_t coding_sequence_size) noexcept {
    return coding_sequence_size % CODON_SIZE;
  }

  /// Construct from a coding sequence and codon index (invalid index yields codon 0, as the reference).
  Codon(const DNA5SequenceCoding& coding_sequence, ContigOffset_t codon_index);

  /// Checked factory: std::nullopt if the codon index is out of range.
  [[nodiscard]] static std::optional<Codon> at(const DNA5SequenceCoding& coding_sequence,
                                               ContigOffset_t codon_index) noexcept;

  /// Construct directly from three symbols.
  [[nodiscard]] static Codon fromSpan(std::span<const Nucleotide, CODON_SIZE> bases) noexcept;

  /// Random access to the codon bases (bounds checked).
  [[nodiscard]] Nucleotide operator[](std::size_t index) const { return bases_.at(index); }

  /// ReturnType the codon as a 3 character base string.
  [[nodiscard]] std::string getSequenceAsString() const;

  /// True if the codon contains the unknown base 'N'.
  [[nodiscard]] bool containsBaseN() const noexcept {
    return bases_[0] == Nucleotide::N or bases_[1] == Nucleotide::N or bases_[2] == Nucleotide::N;
  }

private:

  std::array<Nucleotide, CODON_SIZE> bases_{};

};


}   // end namespace


#endif //KGL_SEQUENCE_CODON_H
