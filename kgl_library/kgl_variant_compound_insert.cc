//
// Created by kellerberrin on 24/11/17.
//


#include "kgl_variant_compound.h"


namespace kgl = kellerberrin::genome;


std::string kgl::CompoundInsert::mutation(char /*delimiter*/, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << "+(" << size() << ")";
  ss << offsetOutput(contigOffset(), output_index);

  return ss.str();

}


bool kgl::CompoundInsert::mutateSequence(SignedOffset_t offset_adjust,
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

  ContigOffset_t remaining_size = dna_sequence_ptr->length() - sequence_offset;
  ContigSize_t reference_size;

  // Check that we are not deleting beyond the end of the sequence.
  if (size() > remaining_size) {

    reference_size = remaining_size;
    ExecEnv::log().vinfo("mutateSequence(), compound insert size: {},  offset: {}, sequence size: {}, remaining size: {}",
                        size(), sequence_offset, dna_sequence_ptr->length(), remaining_size);

  } else {

    reference_size = size();

  }

  auto reference_offset = sequence_offset;

  StringDNA5 seq_string;

  ContigSize_t reference_index = 0;
  for (auto variant : getMap()) {

    std::shared_ptr<const InsertVariant> insert_ptr = std::dynamic_pointer_cast<const InsertVariant>(variant.second);

    if (not insert_ptr) {

      ExecEnv::log().error("mutateSequence(), compound insert contains unexpected variant: {}",
                           variant.second->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }


    // If sufficient sequence remains then check reference.
    if (reference_index < reference_size) {

      // Check the reference.
      if (insert_ptr->reference() != dna_sequence_ptr->at(reference_offset)) {

        ExecEnv::log().info(
        "mutateSequence(), Compound Insert reference base: {} does not match sequence base: {}; Genome: {} Contig: {} Offset: {}",
        DNA5::convertToChar(insert_ptr->reference()),
        DNA5::convertToChar(dna_sequence_ptr->at(reference_offset)),
        insert_ptr->genomeId(), insert_ptr->contig()->contigId(), insert_ptr->offset());

      }

    }


    seq_string.push_back(insert_ptr->mutant());

    ++reference_offset;
    ++reference_index;

  }

  // Mutate the sequence
  DNA5SequenceLinear insert_seq(seq_string);

  if (dna_sequence_ptr->insertSubSequence(sequence_offset, insert_seq)) {

    sequence_size_modify = insert_seq.length();

  } else {

    sequence_size_modify = 0;

  }

  return true;

}
