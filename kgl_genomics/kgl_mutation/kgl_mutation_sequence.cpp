//
// Created by kellerberrin on 7/09/23.
//

#include "kgl_mutation_sequence.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These helper classes keep track of the offset between the zero-offset modified sequence and the contig based offset.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add an indel offset adjust to the map.
// Note that the actual indel increment is passed as an argument.
// But, importantly, the cumulative indel offset is stored in the map.
bool kgl::ModifiedOffsetMap::addModifiedOffset(AdjustedModifiedOffset modified_offset) {

  // Store the cumulative offset;
  SignedOffset_t prev_cumulative_offset = getIndelOffset(modified_offset.contigOffset());
  modified_offset.addCumulativeOffset(prev_cumulative_offset);

  auto [insert_iter, result] = adjust_offset_map_.try_emplace(modified_offset.contigOffset(), modified_offset);
  if (not result) {

    ExecEnv::log().warn("ModifiedOffsetMap::addModifiedOffset; could not insert duplicate offset: {}", modified_offset.contigOffset());

  }

  return result;

}

// Given a contig based sequence interval, return the equivalent zero-based indel modified interval.
kel::OpenRightUnsigned kgl::ModifiedOffsetMap::convertContigModified(const OpenRightUnsigned& contig_interval) const {

  auto [lower, lower_result] = translateZeroOffset(contig_interval.lower());
  auto [upper, upper_result] = translateZeroOffset(contig_interval.upper());

  if (lower_result and upper_result) {

    return {lower, upper};

  }

  ExecEnv::log().warn("ModifiedOffsetMap::convertContigModified; problem converting contig interval: {}", contig_interval.toString());
  return {0, 0};

}

std::pair<kgl::ContigOffset_t, bool> kgl::ModifiedOffsetMap::translateZeroOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the modified sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset) + (getCumulativeOffset() - static_cast<SignedOffset_t>(contig_interval_.lower()));
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("ModifiedOffsetMap::translateZeroOffset; calculated offset: {} is less than zero; contig offset: {}, cumulative indel adjust: {}, contig interval: {}",
                        adjusted_offset, offset, getCumulativeOffset(), contig_interval_.toString());
    return {0, false};

  }

  return {adjusted_offset, true};

}

// Finds the last cumulative offset or 0 if empty map.
kgl::SignedOffset_t kgl::ModifiedOffsetMap::getCumulativeOffset() const {

  if (not adjust_offset_map_.empty()) {

    auto back_iter = adjust_offset_map_.rbegin();
    auto const& [offset, offset_record] = *back_iter;
    return offset_record.indelAdjust();

  }

  return 0;

}

// Finds the map offset element that is <= contig_offset. If not such element exists then returns 0.
kgl::SignedOffset_t kgl::ModifiedOffsetMap::getIndelOffset(ContigOffset_t contig_offset) const {

  auto map_iter = adjust_offset_map_.lower_bound(contig_offset);
  if (map_iter == adjust_offset_map_.end()) {

    if (not adjust_offset_map_.empty()) {

      auto back_iter = adjust_offset_map_.rbegin();
      auto const& [offset, offset_record] = *back_iter;
      return offset_record.indelAdjust();

    }

    return 0;

  }

  auto const& [offset, offset_record] = *map_iter;
  if (offset == contig_offset) {

    return offset_record.indelAdjust();

  }

  // Decrement the iterator to begin() as a guard.
  map_iter = std::ranges::prev(map_iter, 1, adjust_offset_map_.begin());
  auto const& [prev_offset, prev_offset_record] = *map_iter;
  if (prev_offset <= contig_offset) {

    return prev_offset_record.indelAdjust();

  }

  return 0;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class actually modifies the zero-based sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

  auto const [sequence_offset, sequence_result] = translateZeroOffset(interval_update.variantPtr()->offset());
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceSNP; Offset: {} not in interval: {}/{}, update: {}",
                        interval_update.variantPtr()->offset(),
                        contig_interval_.toString(),
                        modified_sequence_.interval().toString(),
                        interval_update.toString());

    return {false, 0};

  }

  const bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, interval_update.variantPtr()->reference());
  if (not reference_match) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceSNP; reference does not match sequence, adjusted offset: {} update: {}",
                        sequence_offset, interval_update.toString());

    return {false, 0};

  }

  if (not modified_sequence_.modifyBase(sequence_offset, interval_update.variantPtr()->alternate().at(0))) {

    ExecEnv::log().error("AdjustedSequence::updateSequenceSNP; could not modify base (SNP) at adjusted offset: {}", sequence_offset);
    return {false, 0};

  }

  return {true, 0};

}

std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceDelete(const SequenceVariantUpdate& interval_update) {

  auto delete_interval = interval_update.priorUpdateIntersect(); // intersection with the prior interval.

  ContigOffset_t ref_offset{0};
  ContigSize_t ref_size{0};

  // Adjust the reference sub-sequence for different deletion types.
  switch(interval_update.updateResult()) {

    case SequenceUpdateResult::PARTIAL_HIGH_DELETE:
      ref_offset = interval_update.variantPtr()->alternateSize();
      ref_size = delete_interval.size();
      break;

    case SequenceUpdateResult::DELETED_REGION:
      ref_offset = interval_update.priorInterval().lower() - interval_update.updatingInterval().lower();
      ref_offset += interval_update.variantPtr()->alternateSize();
      ref_size = interval_update.priorInterval().size();
      break;

    case SequenceUpdateResult::PARTIAL_LOW_DELETE:
      ref_offset = interval_update.updatingInterval().size() - delete_interval.size();
      ref_offset += interval_update.variantPtr()->alternateSize();
      ref_size = interval_update.variantPtr()->referenceSize() - ref_offset;
      break;

    case SequenceUpdateResult::NORMAL:
      ref_offset = interval_update.variantPtr()->alternateSize();
      ref_size = interval_update.variantPtr()->referenceSize() - ref_offset;
      break;

    case SequenceUpdateResult::ERROR:
      ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; encountered bad update: {}", interval_update.toString());
      return {false, 0};

  }

  // Check if the delete extends beyond the interval of interest.
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
  ContigOffset_t contig_delete_offset = interval_update.variantPtr()->offset() + ref_offset;
  auto [sequence_offset, sequence_result] = translateZeroOffset(contig_delete_offset);
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; Offset: {} + reference offset: {}, not in interval: {}/{}, update: {}",
                        interval_update.variantPtr()->offset(),
                        ref_offset,
                        contig_interval_.toString(),
                        modified_sequence_.interval().toString(),
                        interval_update.toString());

    return {false, 0};

  }

  // Explicitly handle the case where the entire region is deleted.
  if (interval_update.updateResult() == SequenceUpdateResult::DELETED_REGION) {

    sequence_offset = 0;

  }

  // Check if the deleted nucleotides match the nucleotides in the modified interval.
  bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, truncated_reference);
  if (not reference_match) {

    // The deleted nucleotides in the modified interval may have already been modified by a prior SNP.
    // So we check the deleted nucleotides against the original unmodified sequence.
    auto const [original_offset, original_result] = originalZeroOffset(interval_update.variantPtr()->offset() + ref_offset);
    reference_match = original_sequence_.compareSubSequence(original_offset, truncated_reference);
    if (not reference_match) {

      OpenRightUnsigned ref_window(sequence_offset, sequence_offset + ref_interval.size());
      const DNA5SequenceLinear interval_reference = modified_sequence_.subSequence(ref_window);

      ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; reference: {} does not match sequence: {} at offset: {}/{}, ref interval: {}, delete interval: {}, Update: {}",
                          truncated_reference.getSequenceAsString(),
                          interval_reference.getSequenceAsString(),
                          sequence_offset,
                          original_offset,
                          ref_interval.toString(),
                          delete_interval.toString(),
                          interval_update.toString()  );

      return {false, 0};

    }

  }

  const OpenRightUnsigned delete_modified{sequence_offset, sequence_offset + ref_size };
  int64_t delete_adjust = -1 * static_cast<int64_t>(ref_size); // deleted size as a -ve.
  bool delete_result = modified_sequence_.deleteSubSequence(delete_modified);
  if (not delete_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; failed to delete : {} from interval: {}",
                        delete_modified.toString(), modified_sequence_.interval().toString());

    return {false, 0};

  }

  bool add_result = modified_offset_map_.addModifiedOffset(AdjustedModifiedOffset(contig_delete_offset, sequence_offset, delete_adjust));
  if (not add_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; failed to add delete adjust: {} at contig offset: {}",
                        delete_adjust, contig_delete_offset);

  }

  return { true, delete_adjust};

}


std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceInsert(const SequenceVariantUpdate& interval_update) {

  ContigOffset_t contig_insert = interval_update.variantPtr()->offset();
  auto const [sequence_offset, sequence_result] = translateZeroOffset(contig_insert);
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; Offset: {}, original: {}/{} not in interval: {}, update: {}",
                        interval_update.variantPtr()->offset(),
                        contig_interval_.toString(),
                        original_sequence_.interval().toString(),
                        modified_sequence_.interval().toString(),
                        interval_update.toString());

    return {false, 0};

  }

  bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, interval_update.variantPtr()->reference());
  if (not reference_match) {

    // The reference letter may have already been modified by a prior SNP.
    // So we check the reference() against the original unmodified sequence.
    auto const [original_offset, original_result] = originalZeroOffset(interval_update.variantPtr()->offset());
    reference_match = original_sequence_.compareSubSequence(original_offset, interval_update.variantPtr()->reference());

    if (not reference_match) {

      ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; reference does not match sequence, adjusted offset: {} update: {}",
                          sequence_offset, interval_update.toString());
      return {false, 0};

    }

  }

  const ContigOffset_t insert_offset = interval_update.variantPtr()->referenceSize();
  contig_insert += insert_offset;  // Adjust the insert offset for the reference size.
  const ContigSize_t insert_size = interval_update.variantPtr()->alternateSize() - insert_offset;

  // Check the insert size is the same size as the updating interval size.
  if (insert_size != interval_update.updatingInterval().size()) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; mismatch between inserted sequence size: {} and update: {}",
                        insert_size, interval_update.toString());

    return {false, 0};

  }

  ContigOffset_t sequence_insert = sequence_offset + insert_offset;
  const DNA5SequenceLinear truncated_alternate = interval_update.variantPtr()->alternate().subSequence(insert_offset, insert_size);

  bool insert_result = modified_sequence_.insertSubSequence(sequence_insert, truncated_alternate);
  if (not insert_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; cannot insert: {} into interval: {}",
                         truncated_alternate.interval().toString(), modified_sequence_.interval().toString());
    return { false, 0};

  }

  bool add_result = modified_offset_map_.addModifiedOffset(AdjustedModifiedOffset(contig_insert,
                                                                                  sequence_offset,
                                                                                  static_cast<SignedOffset_t>(insert_size)));
  if (not add_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; failed to add insert adjust: {} at contig offset: {}",
                        insert_size, sequence_insert);

  }

  return {insert_result, insert_size};

}

// Given an interval map offset, returns the equivalent offset into the zero-offset modified sequence [0, n).
std::pair<kgl::ContigOffset_t, bool> kgl::AdjustedSequence::translateZeroOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the modified sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset) + (modify_offset_adjust_ - static_cast<SignedOffset_t>(contigInterval().lower()));
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("AdjustedSequence::translateZeroOffset; original offset: {}, offset interval: {}, modify_adjust: {}, contig interval: {}, is not in the modified interval: {}",
                        offset, adjusted_offset, modify_offset_adjust_, contigInterval().toString(), modified_sequence_.interval().toString());
    return {0, false};

  }

  if (not modified_sequence_.interval().containsOffset(adjusted_offset)) {

    ExecEnv::log().warn("AdjustedSequence::translateZeroOffset; original offset: {}, offset interval: {} is not in the modified interval: {}",
                        offset, adjusted_offset, modified_sequence_.interval().toString());
    return {0, false};

  }

  return {adjusted_offset, true};

}

// Given an interval map offset, returns the equivalent offset into the zero-offset original sequence [0, m), where (n-m) = base_offset_adjust_.
std::pair<kgl::ContigOffset_t, bool> kgl::AdjustedSequence::originalZeroOffset(ContigOffset_t offset) const {

  // Rebase the offset to the zero offset of the orginal sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(offset)  - static_cast<SignedOffset_t>(contigInterval().lower());
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("AdjustedSequence::originalZeroOffset; offset interval: {} is not in the original interval: {}",
                        adjusted_offset, original_sequence_.interval().toString());
    return {0, false};

  }

  if (not original_sequence_.interval().containsOffset(adjusted_offset)) {

    ExecEnv::log().warn("AdjustedSequence::originalZeroOffset; offset: {} is not in the original interval: {}",
                        adjusted_offset, original_sequence_.interval().toString());
    return {0, false};

  }

  return {adjusted_offset, true};

}
