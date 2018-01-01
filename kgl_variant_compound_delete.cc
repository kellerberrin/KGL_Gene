//
// Created by kellerberrin on 24/11/17.
//



#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;


std::string kgl::CompoundDelete::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;
  ss << "-(" << size() << ")";
  return ss.str() + location(delimiter, output_index);

}


bool kgl::CompoundDelete::mutateCodingSequence(const FeatureIdent_t& sequence_id,
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
  if (not contig()->sequence().offsetWithinSequence(coding_sequence, contigOffset(), sequence_offset,
                                                    sequence_length)) {

    ExecEnv::log().error("mutateCodingSequence(), unable to retrieve coding sequence offset for variant: {}",
                         output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

  // Check that the sequence lengths match.
  if (mutated_sequence->length() != sequence_size) {

    ExecEnv::log().error("compoundDelete::mutateCodingSequence(), sequence length: {} not equal to expected length: {}",
                         mutated_sequence->length(), sequence_size);
    return false;

  }

  SignedOffset_t check_offset = static_cast<SignedOffset_t>(sequence_offset) + offset_adjust;

  std::cout << " adjusted offset:" << check_offset << " seq offset:" << sequence_offset << " adjust:" <<  offset_adjust << std::endl;

  // Check the adjusted offset
  if (check_offset < 0 or check_offset >= static_cast<SignedOffset_t>(mutated_sequence->length())) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {} out of range for sequence length: {}",
                         check_offset, mutated_sequence->length());
    return false;

  }

  ContigOffset_t adjusted_offset = static_cast<ContigOffset_t>(check_offset);

  // check the adjusted sequence offset + the delete size does not exceed the length of the sequence.
  if (adjusted_offset + size() > mutated_sequence->length()) {

    ExecEnv::log().error("mutateCodingSequence(), adjusted offset: {}, with delete size: {} out of range for sequence length: {}",
                         check_offset, size(), mutated_sequence->length());
    return false;

  }

  ContigOffset_t check_reference_offset = adjusted_offset;

  // For all subordinate deletes check that the sequence base code matches the original strand adjusted base code recorded in the variant.
  for (const auto& sub_variant : getMap()) {

    if (mutated_sequence->at(check_reference_offset) != sub_variant.second->strandReference()) {

      ExecEnv::log().warn(
      "mutateCodingSequence(), unexpected; base: {} at seq. offset: {} (strand) reference: {}, probable duplicate variant",
      CodingDNA5::convertToChar(mutated_sequence->at(check_reference_offset)), sequence_offset,
      CodingDNA5::convertToChar(sub_variant.second->strandReference()));

    }

    check_reference_offset++;

  }

  sequence_size_adjust =  -1 * size(); // reduce the sequence size.

  // All is good, so mutate the sequence.
  return mutated_sequence->deleteSubSequence(sequence_offset, size());

}

