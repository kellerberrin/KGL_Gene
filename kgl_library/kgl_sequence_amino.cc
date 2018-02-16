//
// Created by kellerberrin on 31/10/17.
//

#include <set>
#include "kgl_sequence_amino.h"
#include "kgl_genome_db.h"


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



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranslateToAmino - Convert DNA/RNA base sequences to Amino acid sequences.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::shared_ptr<kgl::AminoSequence>
kgl::TranslateToAmino::getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

  StringAminoAcid protein_string;
  AminoAcid::Alphabet amino_acid;

  protein_string.reserve(Codon::codonLength(sequence_ptr));

  for (size_t index = 0; index < Codon::codonLength(sequence_ptr); ++index) {

    amino_acid = table_ptr_->getAmino(Codon(sequence_ptr,index));
    protein_string.push_back(amino_acid);

  }

  std::shared_ptr<AminoSequence> amino_sequence(std::make_shared<AminoSequence>(protein_string));

  return amino_sequence;

}


std::shared_ptr<kgl::AminoSequence>
kgl::TranslateToAmino::getAminoSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                          std::shared_ptr<const DNA5SequenceContig> contig_sequence_ptr) const {

  std::shared_ptr<DNA5SequenceCoding> coding_sequence = contig_sequence_ptr->codingSequence(coding_seq_ptr);
  return getAminoSequence(coding_sequence);

}



size_t kgl::TranslateToAmino::checkNonsenseMutation(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

  for (size_t index = 0; index < Codon::codonLength(sequence_ptr) - 1; ++index) {

    if (table_ptr_->isStopCodon(Codon(sequence_ptr,index))) return index;

  }

  return 0;

}


bool kgl::TranslateToAmino::checkStartCodon(std::shared_ptr<const AminoSequence> sequence_ptr) const {

  if (sequence_ptr->length() > 0) {

    return table_ptr_->isStartAmino(sequence_ptr->at(0));

  } else {

    return false;

  }

}

bool kgl::TranslateToAmino::checkStopCodon(std::shared_ptr<const AminoSequence> sequence_ptr) const {

  if (sequence_ptr->length() > 0) {

    return table_ptr_->isStopAmino(sequence_ptr->at(sequence_ptr->length()-1));

  } else {

    return false;

  }

}


size_t kgl::TranslateToAmino::checkNonsenseMutation(std::shared_ptr<const AminoSequence> sequence_ptr) const {

  for (size_t index = 0; index < sequence_ptr->length() - 1; ++index) {

    if (table_ptr_->isStopAmino(sequence_ptr->at(index))) return index;

  }

  return 0;

}


kgl::AminoAcid::Alphabet kgl::TranslateToAmino::getAmino(const Codon& codon) const {

  if (codon.containsBaseN()) {

    return AminoAcid::Alphabet::Z;

  }

  return table_ptr_->getAmino(codon);

}



kgl::AminoAcid::Alphabet kgl::TranslateToAmino::getAmino(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr,
                                                         ContigOffset_t codon_index) const {

  Codon codon(sequence_ptr, codon_index);

  if (codon.containsBaseN()) {

    return AminoAcid::Alphabet::Z;

  }

  return table_ptr_->getAmino(codon);

}


