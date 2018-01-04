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


bool kgl::CompoundDelete::mutateSequence(SignedOffset_t offset_adjust,
                                         std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) const {

  SignedOffset_t adjusted_offset = offset() + offset_adjust;

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence_ptr->length())) {

    ExecEnv::log().error("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                         adjusted_offset, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);

  auto reference_offset = sequence_offset;

  for (auto variant : getMap()) {

    std::shared_ptr<DeleteVariant const> delete_ptr = std::dynamic_pointer_cast<const DeleteVariant>(variant.second);

    if (not delete_ptr) {

      ExecEnv::log().error("mutateSequence(), compound delete contains unexpected variant: {}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    // Check the reference.
    if (delete_ptr->reference() != dna_sequence_ptr->at(reference_offset)) {

      ExecEnv::log().warn("insertSequence(), Delete reference base: {} does not match sequence base: {} at sequence offset: {}",
                          DNA5::convertToChar(delete_ptr->reference()),
                          DNA5::convertToChar(dna_sequence_ptr->at(reference_offset)),
                          sequence_offset);

    }

    ++reference_offset;

  }

  // Mutate the sequence
  dna_sequence_ptr->deleteSubSequence(sequence_offset, size());

  return true;

}
