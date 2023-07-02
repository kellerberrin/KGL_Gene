//
// Created by kellerberrin on 3/01/18.
//


#include "kgl_variant_mutation_offset.h"


namespace kgl = kellerberrin::genome;



// Important - This routine is paired with updateIndelAccounting(), which should be called AFTER calling (this routine) adjustIndelOffsets()
// Important - calculate indel offset adjustment with (this routine) adjustIndelOffsets() BEFORE calling updateIndelAccounting().
kgl::SignedOffset_t kgl::VariantMutationOffset::adjustIndelOffsets(ContigOffset_t contig_offset) const {

  SignedOffset_t indel_offset_adjust{0};

  // Lookup all indel offsets that are smaller or the same for the variant offset.
  auto upper_bound = indel_accounting_map_.upper_bound(contig_offset);

  for (auto iter = indel_accounting_map_.begin(); iter != upper_bound; ++iter) {

    auto const& [offset, adjust] = *iter;
    indel_offset_adjust += adjust;

  }

  return indel_offset_adjust;

}

// Important - This routine is paired with adjustIndelOffsets(), which should be called BEFORE calling (this routine) updateIndelAccounting()
// Important - update the indel offset accounting structure AFTER the actual indel offset has been calculated with adjustIndelOffsets().
bool kgl::VariantMutationOffset::updateIndelAccounting(ContigOffset_t allele_offset, SignedOffset_t sequence_size_modify) {

  if (sequence_size_modify != 0) {

    auto const [iter, result] = indel_accounting_map_.try_emplace(allele_offset, sequence_size_modify);
    return result;

  }

  return true;

}


kgl::SignedOffset_t kgl::VariantMutationOffset::totalIndelOffset() const {

  SignedOffset_t indel_offset_adjust = 0;

  for (auto const& [offset, adjust] : indel_accounting_map_) {

    indel_offset_adjust += adjust;

  }

  return indel_offset_adjust;

}