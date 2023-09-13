//
// Created by kellerberrin on 7/09/23.
//

#include "kgl_mutation_sequence.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


void kgl::AdjustedSequence::initializeSequences() {

  if (not interval_modify_map_.empty()) {

    auto [offset, first_sequence_modify] = *interval_modify_map_.begin();

    if (contig_interval_ != first_sequence_modify.priorInterval()) {

      ExecEnv::log().error("AdjustedSequence::initializeSequences; contig interval: {} does not match initial map interval: {}",
                           contig_interval_.toString(), first_sequence_modify.priorInterval().toString());
      // Clear the map for good measure.
      interval_modify_map_.clear();

    }

  }

  modified_sequence_ = contig_ref_ptr_->getSubSequence(contig_interval_);
  original_sequence_ = contig_ref_ptr_->getSubSequence(contig_interval_);

}

bool kgl::AdjustedSequence::updateSequence() {

  bool result{true};
  std::pair<bool, SignedOffset_t> update{true, 0};

  for (auto const& [offset, interval_update] : interval_modify_map_) {

    // Adjust the base offset to the current lower offset of the current interval.
    current_offset_adjust_ = interval_update.priorInterval().lower();

    switch(interval_update.variantPtr()->variantType()) {

      case VariantType::SNP:
        update = updateSequenceSNP(interval_update);
        break;

      case VariantType::INDEL_DELETE:
        update = updateSequenceDelete(interval_update);
        break;

      case VariantType::INDEL_INSERT:
        update = updateSequenceInsert(interval_update);
        break;

    }

    auto const& [ update_result, update_offset] = update;

    // Update indel modify offset adjustment.
    modify_offset_adjust_ +=  update_offset;
    result = result and update_result;

  }

  return result;

}

std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceSNP(const SequenceVariantUpdate& interval_update) {

  const ContigOffset_t sequence_offset = translateZeroOffset(interval_update.variantPtr()->offset());

  const bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, interval_update.variantPtr()->reference());
  if (not reference_match) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceSNP; reference does not match sequence, adjusted offset: {} update: {}",
                        sequence_offset, interval_update.toString());

    return {false, 0};

  }

  OpenRightUnsigned prior_interval(modified_sequence_.interval());
  OpenRightUnsigned insert_update(sequence_offset, sequence_offset);

  if (not modified_sequence_.modifyBase(sequence_offset, interval_update.variantPtr()->alternate().at(0))) {

    ExecEnv::log().error("AdjustedSequence::updateSequenceSNP; could not modify base (SNP) at adjusted offset: {}", sequence_offset);
    return {false, 0};

  }

  OpenRightUnsigned post_interval(modified_sequence_.interval());
  SequenceVariantUpdate insert_record(interval_update.variantPtr(),
                                      prior_interval,
                                      post_interval,
                                      insert_update,
                                      SequenceUpdateResult::NORMAL);
  update_audit_map_.insert({interval_update.variantPtr()->offset(), insert_record});

  return {true, 0};

}

std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceDelete(const SequenceVariantUpdate& interval_update) {

  auto delete_interval = interval_update.priorUpdateIntersect(); // intersection with the prior interval.

  // Difference between the size of whole delete interval and the intersection with the prior interval.
  size_t delete_size_adjust = interval_update.updatingInterval().size() - delete_interval.size();

  if (delete_size_adjust > 0) {

    ExecEnv::log().info("AdjustedSequence::updateSequenceDelete: update: {}, truncated: {}/{}, prior interval: {}/{}, modified interval: {}",
                        interval_update.updatingInterval().toString(),
                        delete_interval.toString(),
                        translateZeroOffset(delete_interval).toString(),
                        interval_update.priorInterval().toString(),
                        translateZeroOffset(interval_update.priorInterval()).toString(),
                        modified_sequence_.interval().toString());

    auto trans_prior = translateZeroOffset(interval_update.priorInterval());
    if (trans_prior.lower() > 0) {

      for (auto const& [offset, audit_record] : update_audit_map_) {

        ExecEnv::log().info("AdjustedSequence::updateSequenceDelete; audit: {}", audit_record.toString());

      }

    }

  }

//  return {true, 0};

  // Use the above to calculate the number of deleted nucleotides from the prior interval.
  ContigOffset_t ref_offset = interval_update.variantPtr()->alternateSize() + delete_size_adjust;
  ContigSize_t ref_size = interval_update.variantPtr()->referenceSize() - ref_offset;
  OpenRightUnsigned ref_interval(ref_offset, ref_offset + ref_size);

  // Check the delete size is the same size as the updating interval size.
  if (ref_interval.size() != delete_interval.size()) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; mismatch between delete sequence size: {} and update: {}",
                        ref_interval.size(), delete_interval.toString());

    return {false, 0};

  }

  // Retrieve the interval deleted nucleotides.
  const DNA5SequenceLinear truncated_reference = interval_update.variantPtr()->reference().subSequence(ref_interval);

  // Convert to an equivalent offset in the modified interval.
  const ContigOffset_t sequence_offset = translateZeroOffset(interval_update.variantPtr()->offset() + ref_offset);
  // Check if the deleted nucleotides match the nucleotides in the modified interval.
  bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, truncated_reference);
  if (not reference_match) {

    // The deleted nucleotides in the modified interval may have already been modified by a prior SNP.
    // So we check the deleted nucleotides against the original unmodified sequence.
//    const ContigOffset_t original_offset = originalZeroOffset(delete_interval.lower());
    const ContigOffset_t original_offset = originalZeroOffset(interval_update.variantPtr()->offset() + ref_offset);
    auto orig_reference_match = original_sequence_.compareSubSequence(original_offset, truncated_reference);
    if (not reference_match) {

      if (orig_reference_match) {

        ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; **** original match at: original offset: {}, interval: {}/{} , contig interval: {}, base: {}, current: {}, modify: {}",
                            original_offset,
                            modified_sequence_.interval().toString(),
                            original_sequence_.interval().toString(),
                            contig_interval_.toString(),
                            base_offset_adjust_,
                            current_offset_adjust_,
                            modify_offset_adjust_);

      }

      auto mod_offset = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(sequence_offset)-10);
      OpenRightUnsigned ref_window(mod_offset, mod_offset + 20 + ref_interval.size());
      const DNA5SequenceLinear interval_reference = modified_sequence_.subSequence(ref_window);

      auto mod_original = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(original_offset)-10);
      OpenRightUnsigned orig_window(mod_original, mod_original + ref_interval.size() + 20);
      const DNA5SequenceLinear interval_original = original_sequence_.subSequence(orig_window);

      ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; reference: {}/{} does not match sequence: {}/{}, offset: {}/{} interval: {}, update: {}/{}",
                          truncated_reference.getSequenceAsString(),
                          ref_interval.toString(),
                          interval_reference.getSequenceAsString(),
                          interval_original.getSequenceAsString(),
                          sequence_offset,
                          original_offset,
                          modified_sequence_.interval().toString(),
                          interval_update.toString(),
                          delete_interval.toString());
      return {false, 0};

    }

  }

//  return {true, 0};

  const OpenRightUnsigned delete_modified{sequence_offset, sequence_offset + ref_size };
  bool delete_result = modified_sequence_.deleteSubSequence(delete_modified);
  if (not delete_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; failed to delete : {} from interval: {}",
                        delete_modified.toString(), modified_sequence_.interval().toString());

    return {false, 0};

  }

  return { true, -1 * ref_size};

}

std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceInsert(const SequenceVariantUpdate& interval_update) {

  const ContigOffset_t sequence_offset = translateZeroOffset(interval_update.variantPtr()->offset());
  bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, interval_update.variantPtr()->reference());

  if (not reference_match) {

    // The reference letter may have already been modified by a prior SNP.
    // So we check the reference() against the original unmodified sequence.
    const ContigOffset_t original_offset = originalZeroOffset(interval_update.variantPtr()->offset());
    reference_match = original_sequence_.compareSubSequence(original_offset, interval_update.variantPtr()->reference());

    if (not reference_match) {

      ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; reference does not match sequence, adjusted offset: {} update: {}",
                          sequence_offset, interval_update.toString());
      return {false, 0};

    }

  }

  const ContigOffset_t insert_offset = interval_update.variantPtr()->referenceSize();
  const ContigSize_t insert_size = interval_update.variantPtr()->alternateSize() - insert_offset;

  // Check the insert size is the same size as the updating interval size.
  if (insert_size != interval_update.updatingInterval().size()) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; mismatch between inserted sequence size: {} and update: {}",
                        insert_size, interval_update.toString());

    return {false, 0};

  }

  ContigOffset_t sequence_insert = sequence_offset + insert_offset;

  OpenRightUnsigned prior_interval = translateZeroOffset(interval_update.priorInterval());
  OpenRightUnsigned insert_update(sequence_insert, sequence_insert + insert_size);

  const DNA5SequenceLinear truncated_alternate = interval_update.variantPtr()->alternate().subSequence(insert_offset, insert_size);
  bool insert_result = modified_sequence_.insertSubSequence(sequence_insert, truncated_alternate);
  if (not insert_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; cannot insert: {} into interval: {}",
                         truncated_alternate.interval().toString(), modified_sequence_.interval().toString());
    return { false, 0};

  }

  OpenRightUnsigned post_interval(modified_sequence_.interval());
  SequenceVariantUpdate insert_record(interval_update.variantPtr(),
                                      prior_interval,
                                      post_interval,
                                      insert_update,
                                      SequenceUpdateResult::NORMAL);
  update_audit_map_.insert({interval_update.variantPtr()->offset(), insert_record});

  return {insert_result, insert_size};

}

// Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
kgl::ContigOffset_t kgl::AdjustedSequence::translateZeroOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the modified sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset) + (modify_offset_adjust_ - static_cast<SignedOffset_t>(base_offset_adjust_));
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("AdjustedSequence::translateZeroOffset; offset interval: {} is not in the modified interval: {}",
                        adjusted_offset, modified_sequence_.interval().toString());
    return 0;

  }

  if (not modified_sequence_.interval().containsOffset(adjusted_offset)) {

    ExecEnv::log().warn("AdjustedSequence::translateZeroOffset; offset interval: {} is not in the modified interval: {}",
                        adjusted_offset, modified_sequence_.interval().toString());
    return 0;

  }

  return adjusted_offset;

}

// Same as above but translates an interval.
kel::OpenRightUnsigned kgl::AdjustedSequence::translateZeroOffset(const OpenRightUnsigned& interval) const {

  const SignedOffset_t adjusted_offset =  modify_offset_adjust_ - static_cast<SignedOffset_t>(base_offset_adjust_);
  return interval.translate(adjusted_offset);

}



// Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
kgl::ContigOffset_t kgl::AdjustedSequence::translateRelativeOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the modified sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset) + (modify_offset_adjust_ - static_cast<SignedOffset_t>(current_offset_adjust_));
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("AdjustedSequence::translateRelativeOffset; offset interval: {} is not in the modified interval: {}",
                        adjusted_offset, modified_sequence_.interval().toString());
    return 0;

  }

  if (not modified_sequence_.interval().containsOffset(adjusted_offset)) {

    ExecEnv::log().warn("AdjustedSequence::translateRelativeOffset; offset interval: {} is not in the modified interval: {}",
                        adjusted_offset, modified_sequence_.interval().toString());
    return 0;

  }

  return adjusted_offset;

}

// Same as above but translates an interval.
kel::OpenRightUnsigned kgl::AdjustedSequence::translateRelativeOffset(const OpenRightUnsigned& interval) const {

  const SignedOffset_t adjusted_offset =  modify_offset_adjust_ - static_cast<SignedOffset_t>(current_offset_adjust_);
  return interval.translate(adjusted_offset);

}


// Given an interval map offset, returns the equivalent offset into the zero-offset original sequence [0, m), where (n-m) = base_offset_adjust_.
kgl::ContigOffset_t kgl::AdjustedSequence::originalZeroOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the orginal sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset)  - static_cast<SignedOffset_t>(base_offset_adjust_);
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("AdjustedSequence::originalZeroOffset; offset interval: {} is not in the original interval: {}",
                        adjusted_offset, original_sequence_.interval().toString());
    return 0;

  }

  if (not original_sequence_.interval().containsOffset(adjusted_offset)) {

    ExecEnv::log().warn("AdjustedSequence::originalZeroOffset; offset: {} is not in the original interval: {}",
                        adjusted_offset, original_sequence_.interval().toString());
    return 0;

  }

  return adjusted_offset;

}

// Same as above but translates an interval.
kel::OpenRightUnsigned kgl::AdjustedSequence::translateOriginalOffset(const OpenRightUnsigned& interval) const {

  const SignedOffset_t adjusted_offset = -1 * static_cast<SignedOffset_t>(base_offset_adjust_);
  return interval.translate(adjusted_offset);

}
