//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_amino.h"

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
                                                             ContigOffset_t codon_index) {

  if (codon_index >= codonLength(sequence_ptr)) {

    ExecEnv::log().error("Invalid codon specified index:{}, for coding sequence length:{} (first codon returned)",
                         codon_index, sequence_ptr->length());
    codon_index = 0;
  }

  AminoAcidTypes::Codon codon;
  codon_index = static_cast<ContigOffset_t>(codon_index * 3);
  codon.bases = sequence_ptr->baseAddress(codon_index);
  return codon;

}


kgl::AminoAcidTypes::AminoType kgl::CodingSequenceDNA5::getAmino(std::shared_ptr<DNA5Sequence> sequence_ptr,
                                                                 ContigOffset_t codon_index) const {

  return table_ptr_->getAmino(getCodon(sequence_ptr, codon_index));

}


bool kgl::CodingSequenceDNA5::codonOffset(const SortedCDS& sorted_cds,
                                          ContigOffset_t contig_offset,
                                          ContigOffset_t& codon_offset,
                                          ContigSize_t& base_in_codon) {

  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if (DNA5Sequence::offsetWithinSequence(sorted_cds, contig_offset, sequence_offset, sequence_length)) {

    codon_offset = static_cast<ContigOffset_t>(sequence_offset / 3);
    base_in_codon = static_cast <ContigOffset_t>(sequence_offset % 3);
    return true;

  } else {

    codon_offset = 0;
    base_in_codon = 0;
    return false;

  }

}

bool kgl::CodingSequenceDNA5::SNPMutation(const SortedCDS& sorted_cds,
                                          const std::shared_ptr<const DNA5Sequence>& contig_sequence_ptr,
                                          ContigOffset_t contig_offset,
                                          typename NucleotideColumn_DNA5::NucleotideType reference_base,
                                          typename NucleotideColumn_DNA5::NucleotideType mutant_base,
                                          ContigOffset_t& codon_offset,
                                          typename AminoAcidTypes::AminoType& reference_amino,
                                          typename AminoAcidTypes::AminoType& mutant_amino) const {

  bool result;

  ContigSize_t base_in_codon;
  result = codonOffset(sorted_cds, contig_offset, codon_offset, base_in_codon);

  if (result) {

    auto sequence_offset = static_cast<ContigOffset_t>(codon_offset * 3);
    ContigSize_t codon_size = 3;
    StrandSense strand = sorted_cds.begin()->second->sequence().strand();
    std::shared_ptr<DNA5Sequence> codon_sequence = contig_sequence_ptr->codingSubSequence(sorted_cds,
                                                                                          sequence_offset,
                                                                                          codon_size);
    if (codon_sequence->length() != codon_size) {

      ExecEnv::log().error("SNPMutation(), expected codon sequence size 3: got size: {}", codon_sequence->length());
      result = false;

    }

    switch(strand) {

      case StrandSense::UNKNOWN:
        ExecEnv::log().error("SNPMutation() CDS strands type not set");
        result = false;
        break;

      case StrandSense::FORWARD:
        break;

      case StrandSense::REVERSE:
        reference_base = NucleotideColumn_DNA5::complementNucleotide(reference_base);
        mutant_base = NucleotideColumn_DNA5::complementNucleotide(mutant_base);
        break;

    }

    // Check the reference base in the codon sequence.
    if (result) {

      if (getCodon(codon_sequence, 0).bases[base_in_codon] != reference_base) {

        ExecEnv::log().error("SNPMutation(), reference base: {} does not match codon base: {}",
                             getCodon(codon_sequence, 0).bases[base_in_codon], reference_base);
        result = false;

      }

    }

    // Get the reference and mutant amino acids
    if (result) {

      ExecEnv::log().info("SNPMutation - reference codon: {}, mutation base in codon: {}",
                          codon_sequence->getSequenceString(), base_in_codon);
      reference_amino = getAmino(codon_sequence, 0);
      codon_sequence->modifyBase(mutant_base, base_in_codon);
      ExecEnv::log().info("SNPMutation - mutant codon: {}", codon_sequence->getSequenceString());
      mutant_amino = getAmino(codon_sequence, 0);

    }

  }

  if (not result) {

    codon_offset = 0;
    reference_amino = table_ptr_->getStopAmino();
    mutant_amino = table_ptr_->getStopAmino();
    return false;

  }

  return result;

}
