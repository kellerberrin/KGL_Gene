//
// Created by kellerberrin on 23/04/23.
//

#include "kgl_variant_filter_db.h"

namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters all offsets where there are exactly 2 identical variants (disregarding phase).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::HomozygousFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  if (offset.getVariantArray().size() != 2) {

    return filtered_offset_ptr;

  }
  if (offset.getVariantArray().front()->HGVS() == offset.getVariantArray().back()->HGVS()) {

    filtered_offset_ptr->addVariant(offset.getVariantArray().front());
    filtered_offset_ptr->addVariant(offset.getVariantArray().back());
    return filtered_offset_ptr;

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Ensure max 2 variants per offset.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::DiploidFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  if (offset.getVariantArray().size() <= 2) {

    for (auto const& variant_ptr : offset.getVariantArray()) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

    return filtered_offset_ptr;

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter unique variants disregarding phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::OffsetDB> kgl::UniqueUnphasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants_;
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    auto variant_hash = variant_ptr->HGVS();
    if (not hashed_variants_.contains(variant_hash)) {

      hashed_variants_.insert(variant_hash);
      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter only unique variants including phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::UniquePhasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants_;
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    auto variant_hash = variant_ptr->HGVS_Phase();
    if (not hashed_variants_.contains(variant_hash)) {

      hashed_variants_.insert(variant_hash);
      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}



