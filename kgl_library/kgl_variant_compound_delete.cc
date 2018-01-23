//
// Created by kellerberrin on 24/11/17.
//



#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;


std::string kgl::CompoundDelete::mutation(char /* delimiter */, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << "-(" << size() << ")";
  ss << offsetOutput(contigOffset(), output_index);

  return ss.str();

}

bool kgl::CompoundDelete::mutateSequence(SignedOffset_t offset_adjust,
                                         std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                         SignedOffset_t& sequence_size_modify) const {

  SignedOffset_t adjusted_offset = offset() + offset_adjust;

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence_ptr->length())) {

    ExecEnv::log().error("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                         adjusted_offset, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = dna_sequence_ptr->length() - sequence_offset;
  ContigSize_t delete_size;

  // Check that we are not deleting beyond the end of the sequence.
  if (size() > max_delete_size) {

    delete_size = max_delete_size;
    ExecEnv::log().vinfo("mutateSequence(), compound deletion size: {},  offset: {}, sequence size: {}, max delete size: {}",
                        size(), sequence_offset, dna_sequence_ptr->length(), max_delete_size);

  } else {

    delete_size = size();

  }

  auto reference_offset = sequence_offset;
  ContigSize_t reference_count = 0;

  for (auto variant : getMap()) {

    if (reference_count >= delete_size) break;

    std::shared_ptr<DeleteVariant const> delete_ptr = std::dynamic_pointer_cast<const DeleteVariant>(variant.second);

    if (not delete_ptr) {

      ExecEnv::log().error("mutateSequence(), compound delete contains unexpected variant: {}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    // Check the reference.
    if (delete_ptr->reference() != dna_sequence_ptr->at(reference_offset)) {

      ExecEnv::log().warn("mutateSequence(), delete reference base: {} does not match sequence base: {} at contig: {} offset: {}",
                          DNA5::convertToChar(delete_ptr->reference()),
                          DNA5::convertToChar(dna_sequence_ptr->at(reference_offset)),
                          delete_ptr->contig()->contigId(), delete_ptr->offset());

    }

    ++reference_offset;
    ++reference_count;

  }

  // Mutate the sequence

  if (dna_sequence_ptr->deleteSubSequence(sequence_offset, delete_size)) {

    sequence_size_modify = -1 * delete_size;

  } else {

    sequence_size_modify = 0;

  }

  return true;

}
