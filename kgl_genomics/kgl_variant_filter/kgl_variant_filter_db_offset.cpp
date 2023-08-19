//
// Created by kellerberrin on 11/08/23.
//

#include "kgl_variant_filter_db_offset.h"


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters all offsets where there are identical homozygous variants (disregarding phase).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::HomozygousFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  // Edge condition.
  if (offset.getVariantArray().size() != 2) {

    return filtered_offset_ptr;

  }
  std::map<std::string, std::vector<std::shared_ptr<const Variant>>> variant_map;
  for (auto const& variant_ptr : offset.getVariantArray()) {

    std::string variant_hash = variant_ptr->HGVS();
    if (variant_map.contains(variant_hash)) {

      auto& [hash, vector] = *variant_map.find(variant_hash);
      vector.push_back(variant_ptr);

    } else {

      std::vector<std::shared_ptr<const Variant>> variant_vector{variant_ptr};
      variant_map.try_emplace(variant_hash, variant_vector);

    }

  }

  for (auto const& [hash, vector] : variant_map) {

    if (vector.size() >= 2) {

      for (auto const& variant_ptr : vector) {

        filtered_offset_ptr->addVariant(variant_ptr);

      }

    }

  }

  return filtered_offset_ptr;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters all offsets to only singleton heterozygous variants.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::HeterozygousFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();

  std::map<std::string, std::vector<std::shared_ptr<const Variant>>> variant_map;
  for (auto const& variant_ptr : offset.getVariantArray()) {

    std::string variant_hash = variant_ptr->HGVS();
    if (variant_map.contains(variant_hash)) {

      auto& [hash, vector] = *variant_map.find(variant_hash);
      vector.push_back(variant_ptr);

    } else {

      std::vector<std::shared_ptr<const Variant>> variant_vector{variant_ptr};
      variant_map.try_emplace(variant_hash, variant_vector);

    }

  }

  for (auto const& [hash, vector] : variant_map) {

    if (vector.size() == 1) {

      filtered_offset_ptr->addVariant(vector.front());

    }

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

