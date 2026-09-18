//
// kgl_sequence_amino.cpp — translation and validation facade (merged into AminoTranslationTable).
//

#include "kgl_genetic_code.h"
#include "kel_exec_env.h"


namespace kgl = kellerberrin::genome;


kgl::AminoSequence kgl::AminoTranslationTable::getAminoSequence(const DNA5SequenceCoding& coding_sequence) const {

  const std::size_t sequence_length = coding_sequence.length();
  const std::size_t mod3_length = Codon::codonRemainder(sequence_length);
  if (mod3_length != 0) {

    ExecEnv::log().warn("Coding sequence length: {} is not mod3, mod3(length): {}, protein sequence length: {}",
                        sequence_length, mod3_length, Codon::codonLength(sequence_length));

  }

  std::vector<Amino> protein;
  protein.reserve(Codon::codonLength(sequence_length));

  const auto bases = coding_sequence.getView().span();
  const std::size_t whole_codons = Codon::codonLength(sequence_length);
  for (std::size_t index = 0; index < whole_codons; ++index) {

    protein.push_back(getAmino(Codon::fromSpan(std::span<const Nucleotide, 3>(bases.data() + 3 * index, 3))));

  }

  return AminoSequence(std::move(protein));

}


bool kgl::AminoTranslationTable::checkStartCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {
    return isStartAmino(amino_sequence[0]);
  }
  return false;

}


bool kgl::AminoTranslationTable::checkStopCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {
    return isStopAmino(amino_sequence[amino_sequence.length() - 1]);
  }
  return false;

}


std::pair<std::size_t, bool> kgl::AminoTranslationTable::firstStopSequenceSize(const AminoSequence& amino_sequence) const {

  for (std::size_t index = 0; index < amino_sequence.length(); ++index) {

    if (isStopAmino(amino_sequence[index])) {
      // If found, return the length of the sequence including the stop codon.
      return {index + 1, true};
    }

  }

  // Not found, return the length of the sequence and set the flag to false.
  return {amino_sequence.length(), false};

}
