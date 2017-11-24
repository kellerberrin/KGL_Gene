//
// Created by kellerberrin on 24/11/17.
//


#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;

std::string kgl::CompoundInsert::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << "Compound_Insert>>>>>\n";
  ss << genomeOutput(delimiter, output_index);
  ss << VARIANT_TYPE << delimiter;
  ss << mutation(delimiter, output_index) << "\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output(delimiter, output_index);
  }
  ss << "<<<<<Compound_Insert\n";
  return ss.str();

}

std::string kgl::CompoundInsert::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigSize_t base_in_codon;
    ContigOffset_t codon_offset;
    contig()->sequence().codonOffset(sequence, contigOffset(), codon_offset, base_in_codon);

    ss << "+" << "(" << variant_map_.size() << ")" << offsetOutput(codon_offset, output_index) << delimiter;

  }


  ss << "+" << "(" << variant_map_.size() << ")" << offsetOutput(contigOffset(), output_index) << delimiter;
  return ss.str();

}


bool kgl::CompoundInsert::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                               std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {

  return true;

}

