//
// Created by kellerberrin on 24/11/17.
//



#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;



std::string kgl::CompoundDelete::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << variant_map_.size() << delimiter;
  ss << mutation(delimiter, output_index) << "\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->suboutput(delimiter, output_index);
  }
  return ss.str();

}


std::string kgl::CompoundDelete::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigSize_t base_in_codon;
    ContigOffset_t codon_offset;
    codonOffset(codon_offset, base_in_codon);

    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR
       << offsetOutput(base_in_codon, output_index) << delimiter;

    if (not getMap().empty()) {

      getMap().rbegin()->second->codonOffset(codon_offset, base_in_codon);

      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR
         << offsetOutput(base_in_codon, output_index) << delimiter;

    } else {

      ExecEnv::log().error("Compound delete in contig: {} offset: {} is empty", contigId(), offset());

    }

  } else {

    ExecEnv::log().error("Compound delete in contig: {} offset: {} is not in a coding sequence", contigId(), offset());

  }

  return ss.str();

}


bool kgl::CompoundDelete::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                               std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {

  return true;

}

