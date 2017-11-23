//
// Created by kellerberrin on 5/11/17.
//

#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


std::string kgl::CompoundSNP::output(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index) << delimiter;
  ss << mutation(delimiter, output_index) << delimiter;
  ss << "Compound_SNP>>>>>\n";
  for (const auto& variant : variant_map_) {
    ss << variant.second->output(delimiter, output_index) << "\n";
  }
  ss << "<<<<<Compound_SNP\n";
  return ss.str();

}

std::string kgl::CompoundSNP::mutation(char delimiter, VariantOutputIndex output_index) const {

  std::stringstream ss;
  ss << "-" << "(" << variant_map_.size() << ")" << offsetOutput(contigOffset(), output_index) << delimiter;
  return ss.str();

}

bool kgl::CompoundSNP::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                            std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {

  ExecEnv::log().warn("mutateCodingSequence() not yet implemented for CompoundSNP");
  return false;

}

