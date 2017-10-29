//
// Created by kellerberrin on 29/10/17.
//


#include "kgl_amino_sequence.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::AminoSequence>
kgl::CodingSequenceDNA5::getAminoSequence(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

  typename AminoSequence::ProteinString protein_string;
  AminoAcidTypes::AminoType amino_acid;

  protein_string.reserve(codonLength(sequence_ptr));

  for (size_t index = 0; index < codonLength(sequence_ptr); ++index) {

    amino_acid = table_ptr_->getAmino(getCodon(sequence_ptr,index));
    protein_string.push_back(amino_acid);

  }

  std::shared_ptr<AminoSequence> amino_sequence(std::make_shared<AminoSequence>(protein_string));

  return amino_sequence;

}


size_t kgl::CodingSequenceDNA5::checkNonsenseMutation(std::shared_ptr<DNA5Sequence> sequence_ptr) const {

  for (size_t index = 0; index < codonLength(sequence_ptr) - 1; ++index) {

    if (table_ptr_->isStopCodon(getCodon(sequence_ptr,index))) return index;

  }

  return 0;

}


kgl::AminoAcidTypes::Codon kgl::CodingSequenceDNA5::getCodon(std::shared_ptr<DNA5Sequence> sequence_ptr,
                                                        size_t index) const {

  if (index >= codonLength(sequence_ptr)) {

    ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                         index, sequence_ptr->length());
    index = 0;
  }

  AminoAcidTypes::Codon codon;
  index = index * 3;
  codon.bases = sequence_ptr->baseAddress(index);
  return codon;

}
