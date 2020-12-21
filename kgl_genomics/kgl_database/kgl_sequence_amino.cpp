//
// Created by kellerberrin on 31/10/17.
//

#include <set>
#include "kgl_sequence/kgl_sequence_amino.h"
#include "kgl_database/kgl_genome_db.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Amino Sequence - A container for Amino Acid (protein) sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




bool kgl::AminoSequence::removeTrailingStop() {

  if (alphabet_string_.empty()) {

    ExecEnv::log().warn("Attempt to remove trailing stop amino acid from empty amino sequence");
    return false;

  }

  if (*alphabet_string_.rbegin() == AminoAcid::AMINO_STOP) {

    alphabet_string_.pop_back();
    return true;

  } else {  // Not a stop codon.

    ExecEnv::log().info("removeTrailingStop(). Final amino: {} is not a stop codon.", AminoAcid::convertToChar(*alphabet_string_.rbegin()));
    return false;

  }

}


// Returns a sub-sequence
kgl::AminoSequence kgl::AminoSequence::subSequence( ContigOffset_t sub_sequence_offset, ContigSize_t sub_sequence_length) const {

  AminoSequence sub_sequence;
  if (not getSubsequence(sub_sequence_offset, sub_sequence_length, sub_sequence)) {

    ExecEnv::log().error("AminoSequence::subSequence; Cannot get sub-sequence offset: {} and sub sequence size: {} from sequence length: {}",
                         sub_sequence_offset, sub_sequence_length, length());

  }

  return sub_sequence;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::AminoSequence kgl::TranslateToAmino::getAminoSequence(const DNA5SequenceCoding& coding_sequence) const {

  StringAminoAcid protein_string;
  AminoAcid::Alphabet amino_acid;

  protein_string.reserve(Codon::codonLength(coding_sequence));

  for (size_t index = 0; index < Codon::codonLength(coding_sequence); ++index) {

    amino_acid = table_ptr_->getAmino(Codon(coding_sequence,index));
    protein_string.push_back(amino_acid);

  }

  return AminoSequence(std::move(protein_string));

}


kgl::AminoSequence kgl::TranslateToAmino::getAminoSequence(const std::shared_ptr<const CodingSequence>& coding_seq_ptr,
                                                           const DNA5SequenceContig& contig_sequence) const {

  return getAminoSequence(contig_sequence.codingSequence(coding_seq_ptr));

}



size_t kgl::TranslateToAmino::checkNonsenseMutation(const DNA5SequenceCoding& coding_sequence) const {

  for (size_t index = 0; index < Codon::codonLength(coding_sequence) - 1; ++index) {

    if (table_ptr_->isStopCodon(Codon(coding_sequence,index))) return index;

  }

  return 0;

}


bool kgl::TranslateToAmino::checkStartCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {

    return table_ptr_->isStartAmino(amino_sequence.at(0));

  } else {

    return false;

  }

}

bool kgl::TranslateToAmino::checkStopCodon(const AminoSequence& amino_sequence) const {

  if (amino_sequence.length() > 0) {

    return table_ptr_->isStopAmino(amino_sequence.at(amino_sequence.length()-1));

  } else {

    return false;

  }

}


size_t kgl::TranslateToAmino::checkNonsenseMutation(const AminoSequence& amino_sequence) const {

  for (size_t index = 0; index < amino_sequence.length() - 1; ++index) {

    if (table_ptr_->isStopAmino(amino_sequence.at(index))) return index;

  }

  return 0;

}


kgl::AminoAcid::Alphabet kgl::TranslateToAmino::getAmino(const Codon& codon) const {

  if (codon.containsBaseN()) {

    return AminoAcid::Alphabet::Z;

  }

  return table_ptr_->getAmino(codon);

}



kgl::AminoAcid::Alphabet kgl::TranslateToAmino::getAmino(const DNA5SequenceCoding& coding_sequence,
                                                         ContigOffset_t codon_index) const {

  Codon codon(coding_sequence, codon_index);

  if (codon.containsBaseN()) {

    return AminoAcid::Alphabet::Z;

  }

  return table_ptr_->getAmino(codon);

}


