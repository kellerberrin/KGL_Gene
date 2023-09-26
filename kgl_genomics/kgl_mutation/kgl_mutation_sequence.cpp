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

  auto modified_sequence_opt = contig_ref_ptr->sequence_ptr()->subOptSequence(contigInterval());
  if (not modified_sequence_opt) {

    ExecEnv::log().error("Rrequested modified sub interval: {} out of bounds for contig: {} interval: {}",
                         contigInterval().toString(),
                         contig_ref_ptr->contigId(),
                         contig_ref_ptr->sequence_ptr()->interval().toString());
    return;

  }
  auto original_sequence_opt = contig_ref_ptr->sequence_ptr()->subOptSequence(contigInterval());
  if (not original_sequence_opt) {

    ExecEnv::log().error("Requested original sub interval: {} out of bounds for contig: {} interval: {}",
                         contigInterval().toString(),
                         contig_ref_ptr->contigId(),
                         contig_ref_ptr->sequence_ptr()->interval().toString());
    return;

  }

  modified_sequence_ = std::move(modified_sequence_opt.value());
  original_sequence_ = std::move(original_sequence_opt.value());

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

// The zero-based modified sub-sequence.
// Note that the specified sub_interval must be contained in the initial contig_interval
std::optional<kgl::DNA5SequenceLinear> kgl::AdjustedSequence::modifiedSubSequence(const OpenRightUnsigned& sub_interval) const {

  if (not valid_modified_sequence_) {

    ExecEnv::log().warn("AdjustedSequence::modifiedSubSequence; no valid modified sequence available");
    return std::nullopt;

  }

  // Check sub-interval bounds.
  if (not contigInterval().containsInterval(sub_interval)) {

    ExecEnv::log().warn("Sub interval: {} is not contained in contig interval: {}",
                        sub_interval.toString(), contigInterval().toString());
    return std::nullopt;
  }

  // Convert to the zero-offset modified interval.
  auto [modified_interval, convert_result] = modified_offset_map_.lookupModifiedInterval(sub_interval);
  if (not convert_result) {

    ExecEnv::log().warn("Sub interval: {} cannot be converted to zero offset interval", sub_interval.toString());
    return std::nullopt;

  }

  // Retrieve the modified sub-sequence.
  auto sub_sequence_opt = modified_sequence_.subOptSequence(modified_interval);
  if (not sub_sequence_opt) {

    ExecEnv::log().warn("Sub interval: {} cannot be extracted from modified  sequence: {}",
                        modified_interval.toString(), modified_sequence_.interval().toString());

  }

  return sub_sequence_opt;

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
  auto sub_sequence_opt = original_sequence_.subOptSequence(modified_interval);
  if (not sub_sequence_opt) {

    ExecEnv::log().warn("Sub interval: {} cannot be extracted from modified sequence: {}",
                        modified_interval.toString(), original_sequence_.interval().toString());

  }

  return sub_sequence_opt;

}

