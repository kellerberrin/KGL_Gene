//
// Created by kellerberrin on 5/11/17.
//

#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


std::string kgl::CompoundSNP::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;

  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ContigOffset_t codon_offset;
    AminoAcid::Alphabet reference_amino;
    AminoAcid::Alphabet mutant_amino;

    codonMutation(codon_offset, reference_amino, mutant_amino);

    ss << AminoAcid::convertToChar(reference_amino) << offsetOutput(codon_offset, output_index)
       << AminoAcid::convertToChar(mutant_amino) << delimiter;

  }

  return ss.str();

}


bool kgl::CompoundSNP::codonMutation( ContigOffset_t& codon_offset,
                                      AminoAcid::Alphabet& reference_amino,
                                      AminoAcid::Alphabet& mutant_amino) const {



  if (codingSequences().empty()) {

    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    codon_offset = 0;
    return false;

  }

  ContigSize_t base_in_codon;
  codonOffset(codon_offset, base_in_codon);

  auto sequence_offset = static_cast<ContigOffset_t>(codon_offset * Codon::CODON_SIZE);

  const std::shared_ptr<const CodingSequence> coding_seq_ptr = codingSequences().getFirst();

  std::shared_ptr<DNA5SequenceCoding> codon_sequence = contig()->sequence().codingSubSequence(coding_seq_ptr,
                                                                                              sequence_offset,
                                                                                              Codon::CODON_SIZE);

  if (codon_sequence->length() != Codon::CODON_SIZE) {

    ExecEnv::log().error("codonMutation(), expected codon sequence size 3: got size: {}", codon_sequence->length());
    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    codon_offset = 0;
    return false;

  }

  Codon codon(codon_sequence, 0);  // Create the codon.
  reference_amino = contig()->getAminoAcid(codon);

  for (auto variant : getMap()) {

    std::shared_ptr<const SNPVariant> SNP_ptr = std::dynamic_pointer_cast<const SNPVariant>(variant.second);

    if (not SNP_ptr) {

      ExecEnv::log().error("NON SNP Variant :{} found in Compound SNP",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
      mutant_amino = AminoAcid::AMINO_UNKNOWN;
      codon_offset = 0;
      return false;
    }

    CodingDNA5::Alphabet strand_mutant = SNP_ptr->strandMutant();
    CodingDNA5::Alphabet strand_reference = SNP_ptr->strandReference();

    ContigOffset_t sub_codon_offset;
    SNP_ptr->codonOffset(sub_codon_offset, base_in_codon);

    if (codon_offset != sub_codon_offset) {

      ExecEnv::log().error("codonMutation(), subordinate SNP variant: {} in different codon from compound SNP",
                           SNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
      reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
      mutant_amino = AminoAcid::AMINO_UNKNOWN;
      codon_offset = 0;
      return false;

    }

    if (strand_reference != codon[base_in_codon]) {

      ExecEnv::log().error("codonMutation(), strand reference: {} does not match codon reference: {} for variant: {}",
                           CodingDNA5::convertToChar(strand_reference),
                           CodingDNA5::convertToChar(codon[base_in_codon]),
                           SNP_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
      reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
      mutant_amino = AminoAcid::AMINO_UNKNOWN;
      codon_offset = 0;
      return false;
    }

    codon.modifyBase(base_in_codon, strand_mutant);

  }

  mutant_amino = contig()->getAminoAcid(codon);

  return true;

}


bool kgl::CompoundSNP::mutateSequence(SignedOffset_t offset_adjust,
                                      std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) const {

  for (auto variant : getMap()) {

    std::shared_ptr<SNPVariant const> snp_ptr = std::dynamic_pointer_cast<const SNPVariant>(variant.second);

    if (not snp_ptr) {

      ExecEnv::log().error("mutateSequence(), compound snp contains unexpected variant: {}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    if (not snp_ptr->mutateSequence(offset_adjust, dna_sequence_ptr)) {

      ExecEnv::log().error("mutateSequence(), problem mutating sequence with compound snp: {}",
                           output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

  }

  return true;

}
