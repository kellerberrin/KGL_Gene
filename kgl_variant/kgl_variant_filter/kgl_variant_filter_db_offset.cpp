//
// Created by kellerberrin on 11/08/23.
//

#include "kgl_variant_filter_db_offset.h"


namespace kgl = kellerberrin::genome;


/// Filter to offsets consisting of exactly two identical (homozygous) variants.
std::unique_ptr<kgl::OffsetDB> kgl::HomozygousFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  // Edge condition.
  if (offset.getVariantArray().size() != 2) {

    return filtered_offset_ptr;

  }

  // Group the variants by HGVS hash (unphased identity). A std::map is used deliberately -
  // the filtered variant order is the map's iteration order.
  std::map<std::string, std::vector<std::shared_ptr<const Variant>>> variant_map;
  for (auto const& variant_ptr : offset.getVariantArray()) {

    variant_map[variant_ptr->HGVS()].push_back(variant_ptr);

  }

  // Only identical variant pairs (homozygous) are returned.
  for (auto const& [hash, vector] : variant_map) {

    if (vector.size() >= 2) {

      for (auto const& variant_ptr : vector) {

        filtered_offset_ptr->addVariant(variant_ptr);

      }

    }

  }

  return filtered_offset_ptr;

}


/// Filter to singleton (heterozygous) variants only.
std::unique_ptr<kgl::OffsetDB> kgl::HeterozygousFilter::applyFilter(const OffsetDB& offset) const {

  std::unique_ptr<OffsetDB> filtered_offset_ptr = std::make_unique<OffsetDB>();

  // Group the variants by HGVS hash (unphased identity). A std::map is used deliberately -
  // the filtered variant order is the map's iteration order.
  std::map<std::string, std::vector<std::shared_ptr<const Variant>>> variant_map;
  for (auto const& variant_ptr : offset.getVariantArray()) {

    variant_map[variant_ptr->HGVS()].push_back(variant_ptr);

  }

  // Only unique (singleton, heterozygous) variants are returned.
  for (auto const& [hash, vector] : variant_map) {

    if (vector.size() == 1) {

      filtered_offset_ptr->addVariant(vector.front());

    }

  }

  return filtered_offset_ptr;

}


/// Ensure max 2 variants per offset.
std::unique_ptr<kgl::OffsetDB> kgl::DiploidFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  if (offset.getVariantArray().size() > 2) {

    return filtered_offset_ptr;

  }

  for (auto const& variant_ptr : offset.getVariantArray()) {

    filtered_offset_ptr->addVariant(variant_ptr);

  }

  return filtered_offset_ptr;

}


/// Filter unique variants disregarding phase.
std::unique_ptr<kgl::OffsetDB> kgl::UniqueUnphasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants;
  hashed_variants.reserve(offset.getVariantArray().size());
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    if (hashed_variants.insert(variant_ptr->HGVS()).second) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


/// Filter only unique variants including phase.
std::unique_ptr<kgl::OffsetDB> kgl::UniquePhasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants;
  hashed_variants.reserve(offset.getVariantArray().size());
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    if (hashed_variants.insert(variant_ptr->HGVS_Phase()).second) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}
