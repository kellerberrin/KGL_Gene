//
// Created by kellerberrin on 2/11/17.
//


#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Generate compound delete variants.
//  Given a list of raw delete variants:
// 1. Find delete SNPs in coding sequences.
// 2. Find mod(size, 3) = 0 contiguous SNP deletions.
// 3. If aligned on a codon boundary delete size/3 Amino Acids.
// 4. Or if not aligned, modify 2 Amino Acids and delete size/3 -1 Amino Acids.
// 5. Assemble these SNP deletions into a compound CodonDelete variant.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::CompoundDelete::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index) << delimiter;
  ss << mutation(delimiter, output_index) << "\n";
  ss << "Compound Delete >>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output(delimiter, output_index);
  }
  ss << "<<<<< Compound Delete\n";
  return ss.str();

}

std::string kgl::CompoundDelete::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

  }
  ss << "-" << "(" << variant_map_.size() << ")" << offsetOutput(contigOffset(), output_index) << delimiter;
  return ss.str();

}

bool kgl::CompoundDelete::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                               std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {

  ExecEnv::log().warn("mutateCodingSequence() not yet implemented for CompoundDelete");
  return false;

}

