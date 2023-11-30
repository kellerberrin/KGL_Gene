//
// Created by kellerberrin on 17/09/23.
//

#include "kel_exec_env.h"
#include "kgl_mutation_translate.h"

#include <ranges>

namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These helper classes keep track of the offset between the zero-offset modified sequence and the contig_ref_ptr based offset.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Add an indel offset adjust to the map.
// Note that the actual indel increment is passed as an argument.
// But, importantly, the cumulative indel offset is stored in the map.
bool kgl::ModifiedOffsetMap::addModifiedOffset(AdjustedModifiedOffset modified_offset) {

  // Store the cumulative offset;
  SignedOffset_t cumulative_offset = getIndelOffset(modified_offset.contigOffset()) + modified_offset.offsetAdjust();
  modified_offset.updateCumulative(cumulative_offset);

  auto [insert_iter, result] = adjust_offset_map_.try_emplace(modified_offset.contigOffset(), modified_offset);
  if (not result) {

    ExecEnv::log().warn("ModifiedOffsetMap::addModifiedOffset; could not insert duplicate offset: {}", modified_offset.contigOffset());

  }

  return result;

}


std::string kgl::AdjustedModifiedOffset::toString() const {

  std::string record_string;

  record_string += "Contig Offset: " + std::to_string(contig_offset_);
  record_string += " Zero Offset: " + std::to_string(modified_offset_);
  record_string += " Indel Adjust: " + std::to_string(offset_adjust_);
  record_string += " Cumulative Adjust: " + std::to_string(cumulative_adjust_);
  record_string += " Update Record: " + sequence_update_.toString();

  return record_string;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::pair<kgl::ContigOffset_t, bool> kgl::ModifiedOffsetMap::modifiedZeroOffset(ContigOffset_t contig_offset) const {

  return calcModifiedOffset(contig_offset, getIndelOffset(contig_offset));

}

std::pair<kgl::ContigOffset_t, bool> kgl::ModifiedOffsetMap::calcModifiedOffset(ContigOffset_t contig_offset, SignedOffset_t indel_adjust) const {

  // Rebase the offset to the zero offset of the modified sequence.
  const SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(contig_offset) + indel_adjust - static_cast<SignedOffset_t>(contig_interval_.lower());
  if (adjusted_offset < 0) {

    ExecEnv::log().warn("ModifiedOffsetMap::calcModifiedOffset; calculated offset: {} is less than zero; offset input: {}, cumulative indel adjust: {}, contig_ref_ptr interval: {}",
                        adjusted_offset, contig_offset, indel_adjust, contig_interval_.toString());

    for (auto const& [map_offset, map_record] :  adjust_offset_map_) {

      ExecEnv::log().info("ModifiedOffsetMap::calcModifiedOffset; map offset: {} map record: {}", map_offset, map_record.toString());

    }

    return {0, false};

  }

  return {adjusted_offset, true};

}

// Finds the map offset element that is <  contig_offset. If no such element exists, then returns 0.
kgl::SignedOffset_t kgl::ModifiedOffsetMap::getIndelOffset(ContigOffset_t contig_offset) const {

  auto previous_indel_opt = getPreviousIndel(contig_offset);
  if (not previous_indel_opt) {

    return 0;

  }

  auto previous_indel = previous_indel_opt.value();
  return previous_indel.cumulativeAdjust();

}


// Finds the offset map element that is <  contig_offset. If no such element exists, then returns std::nullopt.
std::optional<kgl::AdjustedModifiedOffset> kgl::ModifiedOffsetMap::getPreviousIndel(ContigOffset_t contig_offset) const {

  auto lower_bound_iter = adjust_offset_map_.lower_bound(contig_offset);
  // Create a sub-range [0, contig_offset).
  auto lb_sub_range = std::ranges::subrange(adjust_offset_map_.begin(), lower_bound_iter);
  // Reverse the sub-range.
  std::ranges::reverse_view reverse_sub_range{lb_sub_range };
  // Search for the first map offset < than the argument contig_offset.
  for (auto const& [map_offset, map_record]: reverse_sub_range) {

    if (map_offset < contig_offset) {

      return AdjustedModifiedOffset(map_record);

    }

  }

  return std::nullopt;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Given a contig_ref_ptr based sequence interval, return the equivalent zero-based unmodified/reference interval.
std::pair<kel::OpenRightUnsigned, bool> kgl::ModifiedOffsetMap::lookupOriginalInterval(const OpenRightUnsigned& contig_interval) const {

  auto [lower, lower_result] = originalZeroOffset(contig_interval.lower());
  auto [upper, upper_result] = originalZeroOffset(contig_interval.upper());

  if (lower_result and upper_result) {

    return {{lower, upper}, true};

  }

  ExecEnv::log().warn("ModifiedOffsetMap::lookupOriginalInterval; problem converting contig_ref_ptr interval: {}", contig_interval.toString());
  return {{0, 0}, false};

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Given a contig_ref_ptr based sequence interval, return the equivalent zero-based indel modified interval.
std::pair<kel::OpenRightUnsigned, bool> kgl::ModifiedOffsetMap::lookupModifiedInterval(const OpenRightUnsigned& contig_interval) const {

  auto [lower, lower_result] = lookupIndelOffset(contig_interval.lower());
  auto [upper, upper_result] = lookupIndelOffset(contig_interval.upper());

  if (lower_result and upper_result) {

    return {{lower, upper}, true};

  }

  ExecEnv::log().warn("ModifiedOffsetMap::lookupModifiedInterval; problem converting contig_ref_ptr interval: {}", contig_interval.toString());
  return {{0, 0}, false};

}


// Lookup an indel modified zero-offset sequence - generally a gene or similar.
// Note that we must account for contig_ref_ptr offsets that occur in the shadow of a delete.
std::pair<kgl::ContigOffset_t, bool> kgl::ModifiedOffsetMap::lookupIndelOffset(ContigOffset_t contig_offset) const {

  auto previous_indel_opt = getPreviousIndel(contig_offset);
  // If no previous indel
  if (not previous_indel_opt) {

    return originalZeroOffset(contig_offset);

  }

  auto& map_record = previous_indel_opt.value();

  // Does the offset fall within a deleted region?
  if (map_record.sequenceUpdate().variantPtr()->variantType() == VariantType::INDEL_DELETE) {

    // Re-construct the deleted interval
    ContigSize_t delete_size = map_record.sequenceUpdate().updatingInterval().size();
    OpenRightUnsigned deleted_interval(map_record.contigOffset(),
                                       map_record.contigOffset() + delete_size);
    // If the offset is within the deleted region.
    // Then adjust the calculated offset to the base of the delete and return.
    if (deleted_interval.containsOffset(contig_offset)) {

      return { map_record.modifiedOffset(), true};

    }

  }

  // Not in the shadow of a delete so just return the calculate modified offset.
  return calcModifiedOffset(contig_offset, map_record.cumulativeAdjust());

}

