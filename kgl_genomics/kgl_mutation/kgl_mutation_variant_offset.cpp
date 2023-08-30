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

bool kgl::AdjustedSequenceOffset::updateOffsetMap(const std::shared_ptr<const Variant>& variant_ptr) {

  VariantType variant_type = variant_ptr->variantType();

  if (variant_type == VariantType::SNP) {

    return true; // No update.

  }

  size_t unmodified_size = modified_interval_.size();
  ContigOffset_t indel_offset{0};

  if (variant_type == VariantType::INDEL_DELETE) {

    indel_offset = variant_ptr->offset() + variant_ptr->alternateSize(); // Actual offset of delete.
    ContigOffset_t adjusted_delete_offset = indel_offset + totalIndelOffset(); // Adjusted for previous indels.
    size_t delete_size = variant_ptr->referenceSize() - variant_ptr->alternateSize(); // Size of the delete.
    OpenRightInterval delete_interval(adjusted_delete_offset, adjusted_delete_offset + delete_size); // Adjusted delete interval
    // Check for that the delete occurs on the modified interval
    if (modified_interval_.disjoint(delete_interval)) {

      ExecEnv::log().warn("AdjustedSequenceOffset::updateOffsetMap; delete variant: {}, delete interval: [{}, {}) disjoint on modified sequence interval [{}, {})",
                          variant_ptr->HGVS(),
                          delete_interval.lower(),
                          delete_interval.upper(),
                          modified_interval_.lower(),
                          modified_interval_.upper());
      return false;

    }

    modified_interval_ = modified_interval_.deleteInterval(delete_interval);

  } else if (variant_type == VariantType::INDEL_INSERT) {

    indel_offset = variant_ptr->offset() + variant_ptr->referenceSize(); // Actual offset of the insert.
    ContigOffset_t adjusted_insert_offset = indel_offset + totalIndelOffset(); // Adjusted for previous indels.
    size_t insert_size = variant_ptr->alternateSize() - variant_ptr->referenceSize(); // Size of the insert.
    OpenRightInterval insert_interval(adjusted_insert_offset, adjusted_insert_offset + insert_size); // Adjusted insert interval
    if (not modified_interval_.containsOffset(adjusted_insert_offset)) {

      ExecEnv::log().warn("AdjustedSequenceOffset::updateOffsetMap; insert variant: {}, insert offset: {} not in sequence interval [{}, {})",
                          variant_ptr->HGVS(),
                          adjusted_insert_offset,
                          modified_interval_.lower(),
                          modified_interval_.upper());
      return false;

    }

    modified_interval_ = modified_interval_.insertInterval(insert_interval);

  }

  size_t modified_size = modified_interval_.size();
  SignedOffset_t signed_change_size = static_cast<SignedOffset_t>(modified_size) - static_cast<SignedOffset_t>(unmodified_size);
  return updateIndelAccounting(indel_offset, signed_change_size);

}


bool kgl::AdjustedSequenceOffset::reconcileIntervalOffset() const {

  return static_cast<SignedOffset_t>(modified_interval_.size())
         == (totalIndelOffset() + static_cast<SignedOffset_t>(original_interval_.size()));

}


