//
// Created by kellerberrin on 3/01/18.
//


#include "kgl_variant_mutation_offset.h"
#include "kgl_variant_mutation.h"


namespace kgl = kellerberrin::genome;



// Important - This routine is paired with updateIndelAccounting(), which should be called AFTER calling (this routine) adjustIndelOffsets()
// Important - calculate indel offset adjustment with (this routine) adjustIndelOffsets() BEFORE calling updateIndelAccounting().
kgl::SignedOffset_t kgl::VariantMutationOffset::adjustIndelOffsets(ContigOffset_t contig_offset) const {

  SignedOffset_t indel_offset_adjust = 0;

  // Lookup all indel offsets that are smaller or the same for the variant offset.
  auto upper_bound = indel_accounting_map_.upper_bound(contig_offset);

  for (auto it = indel_accounting_map_.begin(); it != upper_bound; ++it) {

    indel_offset_adjust += it->second;

  }

  return indel_offset_adjust;

}

// Important - This routine is paired with adjustIndelOffsets(), which should be called BEFORE calling (this routine) updateIndelAccounting()
// Important - update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().
bool kgl::VariantMutationOffset::updateIndelAccounting(std::shared_ptr<const Variant> variant_ptr,
                                                       SignedOffset_t sequence_size_modify) {

  if (sequence_size_modify != 0) {

    std::pair<ContigOffset_t, SignedOffset_t> insert_pair(variant_ptr->offset(), sequence_size_modify);
    auto result = indel_accounting_map_.insert(insert_pair);

    if (not result.second) {

      ExecEnv::log().error("updateIndelAccounting(), Unable to update indel accounting, duplicate offset variant: {}",
                           variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

  }

  return true;

}


kgl::SignedOffset_t kgl::VariantMutationOffset::totalIndelOffset() const {

  SignedOffset_t indel_offset_adjust = 0;

  for (auto it = indel_accounting_map_.begin(); it != indel_accounting_map_.end(); ++it) {

    indel_offset_adjust += it->second;

  }

  return indel_offset_adjust;

}