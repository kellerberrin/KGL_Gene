//
// Created by kellerberrin on 7/09/23.
//

#include "kgl_mutation_sequence.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This class actually modifies the zero-based sequence.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::AdjustedSequence::clear() {

  interval_modify_map_.clear();
  modified_offset_map_.reinitialize({0, 0});
  modified_sequence_.clear();
  original_sequence_.clear();
  valid_modified_sequence_ = false;

}

// Update the unmodified zero-based sequence into the modified sequence.
bool kgl::AdjustedSequence::updateSequence(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                           const OpenRightUnsigned& contig_interval,
                                           const IntervalModifyMap& interval_modify_map) {

  initializeSequences(contig_ref_ptr, contig_interval, interval_modify_map);
  bool result = updateSequence();
  if (not result) {

    ExecEnv::log().warn("djustedSequence::updateSequence; could not update sequence: {}, contig: {}",
                        contig_interval.toString(), contig_ref_ptr->contigId());
    clear();
    return false;

  }

  if (not verifyUpdatedSequence()) {

    ExecEnv::log().warn("djustedSequence::updateSequence; could not verify sequence: {}, contig: {}",
                        contig_interval.toString(), contig_ref_ptr->contigId());
    clear();
    return false;

  }

  valid_modified_sequence_ = true;
  return true;

}


void kgl::AdjustedSequence::initializeSequences(const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                                                const OpenRightUnsigned& contig_interval,
                                                const IntervalModifyMap& interval_modify_map) {

  clear();
  interval_modify_map_ = interval_modify_map;
  modified_offset_map_.reinitialize(contig_interval);

  if (not interval_modify_map_.empty()) {

    auto [offset, first_sequence_modify] = *interval_modify_map_.begin();

    if (contigInterval() != first_sequence_modify.priorInterval()) {

      ExecEnv::log().error("AdjustedSequence::initializeSequences; contig interval: {} does not match initial map interval: {}",
                           contigInterval().toString(), first_sequence_modify.priorInterval().toString());
      // Clear the map for good measure.
      clear();
      return;

    }

  }

  modified_sequence_ = contig_ref_ptr->getSubSequence(contigInterval());
  original_sequence_ = contig_ref_ptr->getSubSequence(contigInterval());

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

    // Update indel modify offset adjustment.
    auto const& [ update_result, update_offset] = update;
    result = result and update_result;

  }

  return result;

}


bool kgl::AdjustedSequence::verifyUpdatedSequence() const {

  bool result{true};
  for (auto const& [offset, interval_update] : interval_modify_map_) {

    const std::shared_ptr<const Variant>& variant_ptr = interval_update.variantPtr();
    if (variant_ptr->variantType() == VariantType::SNP) {

      auto [modified_offset, offset_result] = sequenceMap().modifiedZeroOffset(offset);
      if (not offset_result) {

        ExecEnv::log().warn("AdjustedSequence::verifyUpdatedSequence; Cannot convert offset: {} to modified offset, variant: {}",
                            offset, variant_ptr->HGVS());
        result = false;
        continue;

      }

      const bool reference_match = modified_sequence_.compareSubSequence(modified_offset, variant_ptr->alternate());
      if (not reference_match) {

        ExecEnv::log().warn("AdjustedSequence::verifyUpdatedSequence; alternate does not match modified sequence, modified offset: {} variant: {}",
                            offset, variant_ptr->HGVS());

        result = false;
        continue;

      }

    }

  }

  return result;

}


std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceSNP(const SequenceVariantUpdate& interval_update) {

  ContigOffset_t offset_SNP = interval_update.variantPtr()->offset();
  auto [sequence_offset, sequence_result] = modified_offset_map_.modifiedZeroOffset(offset_SNP);
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceSNP; Offset: {} not in interval: {}/{}, update: {}",
                        interval_update.variantPtr()->offset(),
                        contigInterval().toString(),
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
  auto [sequence_offset, sequence_result] = modified_offset_map_.modifiedZeroOffset(contig_delete_offset);
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; Offset: {} + reference offset: {}, not in interval: {}/{}, update: {}",
                        interval_update.variantPtr()->offset(),
                        ref_offset,
                        contigInterval().toString(),
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
    auto [original_offset, original_result] = modified_offset_map_.originalZeroOffset(contig_delete_offset);
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

  bool add_result = modified_offset_map_.addModifiedOffset(AdjustedModifiedOffset(contig_delete_offset,
                                                                                  sequence_offset,
                                                                                  delete_adjust,
                                                                                  interval_update));
  if (not add_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceDelete; failed to add delete adjust: {} at contig offset: {}",
                        delete_adjust, contig_delete_offset);

  }

  return { true, delete_adjust};

}


std::pair<bool, kgl::SignedOffset_t> kgl::AdjustedSequence::updateSequenceInsert(const SequenceVariantUpdate& interval_update) {

  ContigOffset_t contig_insert = interval_update.variantPtr()->offset();
  auto [sequence_offset, sequence_result] = modified_offset_map_.modifiedZeroOffset(contig_insert);
  if (not sequence_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; Offset: {}, original: {}/{} not in interval: {}, update: {}",
                        interval_update.variantPtr()->offset(),
                        contigInterval().toString(),
                        original_sequence_.interval().toString(),
                        modified_sequence_.interval().toString(),
                        interval_update.toString());

    return {false, 0};

  }

  bool reference_match = modified_sequence_.compareSubSequence(sequence_offset, interval_update.variantPtr()->reference());
  if (not reference_match) {

    // The reference letter may have already been modified by a prior SNP.
    // So we check the reference() against the original unmodified sequence.
    auto [original_offset, original_result] = modified_offset_map_.originalZeroOffset(contig_insert);
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
                                                                                  static_cast<SignedOffset_t>(insert_size),
                                                                                  interval_update));
  if (not add_result) {

    ExecEnv::log().warn("AdjustedSequence::updateSequenceInsert; failed to add insert adjust: {} at contig offset: {}",
                        insert_size, sequence_insert);

  }

  return {insert_result, insert_size};

}

// The zero-based modified sub-sequence.
// Note that the specified sub_interval must be contained in the initial contig_interval
std::optional<kgl::DNA5SequenceLinear> kgl::AdjustedSequence::modifiedSubSequence(const OpenRightUnsigned& sub_interval) const {

  if (not valid_modified_sequence_) {

    ExecEnv::log().warn("AdjustedSequence::modifiedSubSequence; no valid modified sequence available");
    return std::nullopt;

  }

  // Check sub-interval bounds.
  if (not contigInterval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("AdjustedSequence::modifiedSubSequence; sub interval: {} is not contained in contig interval: {}",
                        sub_interval.toString(), contigInterval().toString());
    return std::nullopt;
  }

  // Convert to the zero-offset modified interval.
  auto [modified_interval, convert_result] = modified_offset_map_.lookupModifiedInterval(sub_interval);
  if (not convert_result) {

    ExecEnv::log().warn("AdjustedSequence::modifiedSubSequence; sub interval: {} cannot be converted to zero offset interval",
                        sub_interval.toString());
    return std::nullopt;

  }

  // Retrieve the modified sub-sequence.
  DNA5SequenceLinear sub_sequence = modified_sequence_.subSequence(modified_interval);

  return sub_sequence;

}

// The zero-based unmodified sub-sequence.
// Note that the specified sub_interval must be contained in the initial contig_interval
std::optional<kgl::DNA5SequenceLinear> kgl::AdjustedSequence::originalSubSequence(const OpenRightUnsigned& sub_interval) const {

  if (not valid_modified_sequence_) {

    ExecEnv::log().warn("AdjustedSequence::modifiedSubSequence; no valid original sequence available");
    return std::nullopt;

  }

  // Check sub-interval bounds.
  if (not contigInterval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("AdjustedSequence::originalSubSequence; sub interval: {} is not contained in contig interval: {}",
                        sub_interval.toString(), contigInterval().toString());
    return std::nullopt;
  }

  // Convert to the zero-offset original interval.
  auto [modified_interval, convert_result] = modified_offset_map_.lookupOriginalInterval(sub_interval);
  if (not convert_result) {

    ExecEnv::log().warn("AdjustedSequence::originalSubSequence; sub interval: {} cannot be converted to zero offset interval",
                        sub_interval.toString());
    return std::nullopt;

  }

  // Retrieve the modified sub-sequence.
  DNA5SequenceLinear sub_sequence = original_sequence_.subSequence(modified_interval);

  return sub_sequence;

}

