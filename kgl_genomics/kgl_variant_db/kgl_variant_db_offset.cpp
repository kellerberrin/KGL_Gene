//
// Created by kellerberrin on 4/2/21.
//


#include "kgl_variant_db_offset.h"
#include "kgl_variant_filter_db_genome.h"


namespace kgl = kellerberrin::genome;


std::unique_ptr<kgl::OffsetDB> kgl::OffsetDB::viewFilter(const BaseFilter &filter) const {


  if (filter.filterType() == FilterBaseType::OFFSET_FILTER) {

    const FilterOffsets& offset_filter = static_cast<const FilterOffsets&>(filter);
    return offset_filter.applyFilter(*this);

  }


  // Filter the variants.
  std::unique_ptr<OffsetDB> filtered_offset_ptr(std::make_unique<OffsetDB>());

  // Only variant filters should be at this level, so we check for this.
  if (filter.filterType() != FilterBaseType::VARIANT_FILTER) {

    ExecEnv::log().error("OffsetDB::viewFilter; Filter: {} is not a variant_filter.", filter.filterName());
    return filtered_offset_ptr;

  }

  // Hoist the filter cast out of the variant loop. Note that variant filters must be stateless,
  // the filter object is applied concurrently to all variants without cloning.
  const auto& variant_filter = static_cast<const FilterVariants&>(filter);

  for (const auto& variant_ptr : variant_vector_) {

    if (variant_filter.applyFilter(*variant_ptr)) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


std::pair<size_t, size_t> kgl::OffsetDB::selfFilter(const BaseFilter &filter) {

  // Offset filters are applied to the whole offset. The filter is cloned before use because
  // filters may maintain (mutable) state that must not transfer between filter invocations.
  if (filter.filterType() == FilterBaseType::OFFSET_FILTER) {

    std::shared_ptr<const FilterOffsets> offset_filter = std::dynamic_pointer_cast<const FilterOffsets>(filter.clone());
    auto filtered_offset_ptr = offset_filter->applyFilter(*this);

    size_t prior_count = variant_vector_.size();
    variant_vector_ = std::move(filtered_offset_ptr->variant_vector_);
    size_t post_count = variant_vector_.size();

    return { prior_count, post_count};

  }

  // Only variant filters should be at this level, so we check for this.
  if (filter.filterType() != FilterBaseType::VARIANT_FILTER) {

    ExecEnv::log().error("OffsetDB::selfFilter; Filter: {} is not a variant_filter.", filter.filterName());
    return { variant_vector_.size(), variant_vector_.size() };

  }

  // Filter the variants.
  // Hoist the filter cast out of the variant loop. Note that variant filters must be stateless,
  // the filter object is applied concurrently to all variants without cloning.
  const auto& variant_filter = static_cast<const FilterVariants&>(filter);

  OffsetDBArray filtered_variants;
  for (const auto& variant_ptr : variant_vector_) {

    if (variant_filter.applyFilter(*variant_ptr)) {

      filtered_variants.push_back(variant_ptr);

    }

  }

  size_t prior_count = variant_vector_.size();
  variant_vector_ = std::move(filtered_variants);

  return { prior_count, variant_vector_.size() };

}