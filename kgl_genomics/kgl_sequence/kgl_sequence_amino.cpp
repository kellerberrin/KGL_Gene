//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_sequence_amino.h"


namespace kgl = kellerberrin::genome;




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::AminoSequence kgl::TranslateToAmino::getAminoSequence(const DNA5SequenceCoding& coding_sequence) const {

  StringAminoAcid protein_string;

  size_t sequence_length = coding_sequence.length();
  size_t mod3_length = Codon::codonRemainder(sequence_length);
  if (mod3_length != 0) {

    ExecEnv::log().warn("Coding sequence length: {} is not mod3, mod3(length): {}, protein sequence length: {}",
                        sequence_length, mod3_length, Codon::codonLength(coding_sequence));

  }

  protein_string.reserve(Codon::codonLength(coding_sequence));

  for (size_t index = 0; index < Codon::codonLength(coding_sequence); ++index) {

    protein_string.push_back(getAmino(Codon(coding_sequence,index)));

  }

  return AminoSequence(std::move(protein_string));

}


bool kgl::TranslateToAmino::checkStartCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {

    return table_.isStartAmino(amino_sequence.at(0));

  } else {

    return false;

  }

}

bool kgl::TranslateToAmino::checkStopCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {

    return table_.isStopAmino(amino_sequence.at(amino_sequence.length()-1));

  } else {

    return false;

  }

}


std::pair<size_t, bool> kgl::TranslateToAmino::firstStopSequenceSize(const AminoSequence& amino_sequence) const {

  for (size_t index = 0; index < amino_sequence.length(); ++index) {

    if (table_.isStopAmino(amino_sequence.at(index))) {

      // If found, return the length of the sequence including the stop codon.
      return {index + 1, true};

    }

  }

  // Not found, return the length of the sequence and set the flag to false.
  return { amino_sequence.length(), false};

}


kgl::AminoAcid::Alphabet kgl::TranslateToAmino::getAmino(const Codon& codon) const {

  return table_.getAmino(codon);

}