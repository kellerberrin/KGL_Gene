//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_single.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::SNPVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
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



bool kgl::SNPVariant::mutateCodingSequence(const FeatureIdent_t& sequence_id,
                                           SignedOffset_t offset_adjust,  // Adjust the variant offsets before mutation
                                           ContigSize_t sequence_size,  // Calculated sequence size before mutation.
                                           SignedOffset_t& sequence_size_adjust,  // How the variant modifies sequence size.
                                           std::shared_ptr<DNA5SequenceCoding>& mutated_sequence) const {


  CodingSequenceArray coding_sequence_array = codingSequences();

  // Check that we have a variant in a coding sequence.
  if (coding_sequence_array.empty()) {

    ExecEnv::log().warn("mutateCodingSequence(), variant: {} not in a coding sequence",
                        output(' ', VariantOutputIndex::START_0_BASED, true));
    return true; // just ignored.

  }

  std::shared_ptr<const CodingSequence> coding_sequence = coding_sequence_array.getFirst();

  // Check the sequence id.
  if (coding_sequence->getCDSParent()->id() != sequence_id) {

    ExecEnv::log().warn("mutateCodingSequence(), variant: {} does not mutate sequence id: {}",
                        output(' ', VariantOutputIndex::START_0_BASED, true), sequence_id);
    return true; // just ignored.

  }

  // Get the variant sequence offset
  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if(not contig()->sequence().offsetWithinSequence(coding_sequence, contigOffset(), sequence_offset, sequence_length)) {

    ExecEnv::log().error("mutateCodingSequence(), unable to retrieve coding sequence offset for variant: {}",
                         output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

  // Check that the sequence lengths match.
  if (sequence_length != sequence_size) {

    ExecEnv::log().error("mutateCodingSequence(), unexpected; variant sequence length: {} not equal mutate length: {}",
                         sequence_length, mutated_sequence->length());
    return false;

  }

  SignedOffset_t check_offset = static_cast<SignedOffset_t>(sequence_offset) + offset_adjust;

  // Check the adjusted offset
  if (check_offset < 0 or check_offset >= static_cast<SignedOffset_t>(mutated_sequence->length())) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {} out of range for sequence length: {}",
                        check_offset, mutated_sequence->length());
    return false;

  }

  ContigOffset_t adjusted_offset = static_cast<ContigOffset_t>(check_offset);

  // Check that the sequence base code matches the original strand adjusted base code recorded in the variant.
  if (mutated_sequence->at(adjusted_offset) != strandReference()) {

    ExecEnv::log().warn("mutateCodingSequence(), unexpected; base: {} at seq. offset: {} (strand) reference: {}, probable duplicate variant",
    CodingDNA5::convertToChar(mutated_sequence->at(sequence_offset)), sequence_offset,
    CodingDNA5::convertToChar(strandReference()));

  }

  sequence_size_adjust = 0;

  // All is good, so mutate the sequence.
  return mutated_sequence->modifyBase(sequence_offset, strandMutant());


}

