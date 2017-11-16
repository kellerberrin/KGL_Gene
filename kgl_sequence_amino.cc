//
// Created by kellerberrin on 31/10/17.
//

#include <set>
#include "kgl_sequence_amino.h"
#include "kgl_genome_db.h"


namespace kgl = kellerberrin::genome;


kgl::ProteinString kgl::AminoSequence::emphasizeProteinString(const ProteinString& protein_string,
                                                              const std::vector<ContigOffset_t>& emphasize_offsets) {

  if (emphasize_offsets.empty()) return protein_string;

  // Order the offsets. A < B < C
  std::set<ContigOffset_t> ordered_offsets;
  for (auto offset : emphasize_offsets) {

    ordered_offsets.insert(offset);

  }

  ProteinString emph_protein_string;
  size_t index = 0;
  for (auto offset : ordered_offsets) {

    if (offset >= protein_string.length()) {

      ExecEnv::log().error("emphasizeProteinString() emphasize offset: {} >= protein string length: {}",
                           offset, emph_protein_string.length());
      break;
    }

    while (index < offset) {

      emph_protein_string += protein_string[index];
      index++;

    }

    // Add spaces to emphasize.
    if (index == offset) {

      emph_protein_string += ' ';
      emph_protein_string += protein_string[index];;
      emph_protein_string += ' ';
      index++;

    }

  }

  // Add in the rest of the sequence
  while (index < protein_string.length()) {

    emph_protein_string += protein_string[index];
    index++;

  }

  return emph_protein_string;

}


std::shared_ptr<kgl::AminoSequence>
kgl::TranslateToAmino::getAminoSequence(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

  ProteinString protein_string;
  Amino_t amino_acid;

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



size_t kgl::TranslateToAmino::checkNonsenseMutation(std::shared_ptr<DNA5SequenceCoding> sequence_ptr) const {

  for (size_t index = 0; index < Codon::codonLength(sequence_ptr) - 1; ++index) {

    if (table_ptr_->isStopCodon(Codon(sequence_ptr,index))) return index;

  }

  return 0;

}



kgl::Amino_t kgl::TranslateToAmino::getAmino(std::shared_ptr<DNA5SequenceCoding> sequence_ptr,
                                             ContigOffset_t codon_index) const {

  Codon codon(sequence_ptr, codon_index);

  if (codon.containsBaseN()) {

    return AminoAcid::UNKNOWN_AMINO;

  }

  return table_ptr_->getAmino(codon);

}


bool kgl::TranslateToAmino::SNPMutation(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                          const std::shared_ptr<const DNA5SequenceContig>& contig_sequence_ptr,
                                          ContigOffset_t contig_offset,
                                          DNA5::Alphabet reference_base,
                                          DNA5::Alphabet mutant_base,
                                          ContigOffset_t& codon_offset,
                                          Amino_t& reference_amino,
                                          Amino_t& mutant_amino) const {

  bool result;

  ContigSize_t base_in_codon;
  result = contig_sequence_ptr->codonOffset(coding_seq_ptr, contig_offset, codon_offset, base_in_codon);

  if (result) {

    auto sequence_offset = static_cast<ContigOffset_t>(codon_offset * 3);
    ContigSize_t codon_size = 3;
    StrandSense strand = coding_seq_ptr->getCDSParent()->sequence().strand();
    std::shared_ptr<DNA5SequenceCoding> codon_sequence = contig_sequence_ptr->codingSubSequence(coding_seq_ptr,
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
        reference_base = DNA5::complementNucleotide(reference_base);
        mutant_base = DNA5::complementNucleotide(mutant_base);
        break;

    }

    // Check the reference base in the codon sequence.
    if (result) {

      if (Codon(codon_sequence, 0)[base_in_codon] != reference_base) {

        ExecEnv::log().error("SNPMutation(), reference base: {} does not match codon base: {}",
                             DNA5::convertToChar(Codon(codon_sequence, 0)[base_in_codon]),
                             DNA5::convertToChar(reference_base));
        result = false;

      }

    }

    // Get the reference and mutant amino acids
    if (result) {

      reference_amino = getAmino(codon_sequence, 0);
      codon_sequence->modifyBase(mutant_base, base_in_codon);
      mutant_amino = getAmino(codon_sequence, 0);

    }

  } else {



  }

  if (not result) {

    codon_offset = 0;
    reference_amino = table_ptr_->getStopAmino();
    mutant_amino = table_ptr_->getStopAmino();
    return false;

  }

  return result;

}

