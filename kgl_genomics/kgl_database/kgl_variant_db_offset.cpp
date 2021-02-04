//
// Created by kellerberrin on 4/2/21.
//


#include "kgl_variant_db_offset.h"


namespace kgl = kellerberrin::genome;


std::pair<size_t, size_t> kgl::OffsetDB::inSituFilter(const VariantFilter &filter) {

  std::pair<size_t, size_t> offset_count{0, 0};

  offset_count.first += variant_vector_.size();
  OffsetDBArray filtered_variants;
  filtered_variants.reserve(offset_count.first);

  for (auto const &variant_ptr : variant_vector_) {

    if (filter.applyFilter(*variant_ptr)) {

      filtered_variants.push_back(variant_ptr);

    }

  }

  offset_count.second = filtered_variants.size();

  if (offset_count.first != offset_count.second) {

    if (not filtered_variants.empty()) {

      variant_vector_ = std::move(filtered_variants);


    } else {

      variant_vector_.clear();

    }

  }

  return offset_count;

}