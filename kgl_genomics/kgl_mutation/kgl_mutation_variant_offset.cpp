//
// Created by kellerberrin on 3/01/18.
//


#include "kgl_mutation_variant_offset.h"


namespace kgl = kellerberrin::genome;



// Important - This routine is paired with updateIndelAccounting(), which should be called AFTER calling (this routine) adjustIndelOffsets()
// Important - calculate indel offset adjustment with (this routine) adjustIndelOffsets() BEFORE calling updateIndelAccounting().
kgl::SignedOffset_t kgl::AdjustedSequenceOffset::adjustIndelOffsets(ContigOffset_t contig_offset) const {

  SignedOffset_t indel_offset_adjust{0};

  // Lookup all indel offsets that are smaller or the same for the variant offset.
  auto const upper_bound = indel_accounting_map_.upper_bound(contig_offset);
  auto const lower_bound = indel_accounting_map_.begin();

  for (auto const& [offset, adjust] : std::ranges::subrange(lower_bound, upper_bound)) {

    indel_offset_adjust += adjust;

  }

  return indel_offset_adjust;

}

// Important - This routine is paired with adjustIndelOffsets(), which should be called BEFORE calling (this routine) updateIndelAccounting()
// Important - update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().
bool kgl::AdjustedSequenceOffset::updateIndelAccounting(ContigOffset_t allele_offset, SignedOffset_t sequence_size_modify) {

  bool result{true};

  auto const [iter, insert_result] = indel_accounting_map_.try_emplace(allele_offset, sequence_size_modify);
  // Check that the interval has the same size as the accounting map.
  if (not insert_result) {

    ExecEnv::log().warn("AdjustedSequenceOffset::updateIndelAccounting; cannot update (duplicate) accounting map, oiffset: {}, size:{}",
                        allele_offset, sequence_size_modify);
    result = false;

  }

  if (not reconcileIntervalOffset()) {

    ExecEnv::log().warn("AdjustedSequenceOffset::updateIndelAccounting; accounting map size: {} does not equal calculated interval size: {}",
                        totalIndelOffset() , modified_interval_.size());
    result = false;

  }

  return result;

}


kgl::SignedOffset_t kgl::AdjustedSequenceOffset::totalIndelOffset() const {

  SignedOffset_t indel_offset_adjust{0};

  for (auto const& [offset, adjust] : indel_accounting_map_) {

    indel_offset_adjust += adjust;

  }

  return indel_offset_adjust;

}



bool kgl::AdjustedSequenceOffset::reconcileIntervalOffset() const {

  return static_cast<SignedOffset_t>(modified_interval_.size())
         == (totalIndelOffset() + static_cast<SignedOffset_t>(original_interval_.size()));

}


