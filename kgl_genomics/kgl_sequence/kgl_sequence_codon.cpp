//
// kgl_sequence_codon.cpp — Codon construction and formatting.
//

#include "kgl_sequence_codon.h"
#include "kgl_sequence_dna.h"
#include "kel_exec_env.h"

#include <format>


namespace kgl = kellerberrin::genome;


kgl::ContigSize_t kgl::Codon::codonLength(const DNA5SequenceCoding& coding_sequence) noexcept {

  return codonLength(coding_sequence.length());

}


kgl::Codon::Codon(const DNA5SequenceCoding& coding_sequence, ContigOffset_t codon_index) {

  if (codon_index >= codonLength(coding_sequence.length())) {

    ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                         codon_index, coding_sequence.length());
    codon_index = 0;

  }

  const ContigOffset_t base_offset = codon_index * CODON_SIZE;
  bases_[0] = coding_sequence[base_offset];
  bases_[1] = coding_sequence[base_offset + 1];
  bases_[2] = coding_sequence[base_offset + 2];

}


std::optional<kgl::Codon> kgl::Codon::at(const DNA5SequenceCoding& coding_sequence,
                                         ContigOffset_t codon_index) noexcept {

  if (codon_index >= codonLength(coding_sequence.length())) {
    return std::nullopt;
  }

  const ContigOffset_t base_offset = codon_index * CODON_SIZE;
  Codon codon;
  codon.bases_[0] = coding_sequence[base_offset];
  codon.bases_[1] = coding_sequence[base_offset + 1];
  codon.bases_[2] = coding_sequence[base_offset + 2];
  return codon;

}


kgl::Codon kgl::Codon::fromSpan(std::span<const Nucleotide, CODON_SIZE> bases) noexcept {

  Codon codon;
  codon.bases_[0] = bases[0];
  codon.bases_[1] = bases[1];
  codon.bases_[2] = bases[2];
  return codon;

}


std::string kgl::Codon::getSequenceAsString() const {

  return std::format("{}{}{}",
                     CodingDNA5::convertToChar(bases_[0]),
                     CodingDNA5::convertToChar(bases_[1]),
                     CodingDNA5::convertToChar(bases_[2]));

}
