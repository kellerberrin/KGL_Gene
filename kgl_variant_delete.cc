//
// Created by kellerberrin on 23/12/17.
//


#include "kgl_variant_single.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// InsertVariant - generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::DeleteVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << quality() << delimiter;
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence()->output(delimiter, output_index);

  }

  ss << '\n';

  return ss.str();

}


// This mutates a coding sequence that has already been generated using a CodingSequence (CDS) object.
bool kgl::DeleteVariant::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                              SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                                              ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                                              SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                                              std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {


  sequence_size_adjust = 0;
  return true;

}
