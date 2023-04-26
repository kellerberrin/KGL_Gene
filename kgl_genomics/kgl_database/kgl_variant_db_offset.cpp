//
// Created by kellerberrin on 4/2/21.
//


#include "kgl_variant_db_offset.h"
#include "kgl_variant_filter_db.h"
#include "kgl_variant_filter.h"


namespace kgl = kellerberrin::genome;


std::unique_ptr<kgl::OffsetDB> kgl::OffsetDB::copyFilter(const BaseFilter &filter) const {


  std::shared_ptr<const FilterOffsets> offset_filter = std::dynamic_pointer_cast<const FilterOffsets>(filter.clone());
  if (offset_filter) {

    return offset_filter->applyFilter(*this);

  }

  // Only variant filters should be at this level, so we check for this.
  std::shared_ptr<const FilterVariants> variant_filter = std::dynamic_pointer_cast<const FilterVariants>(filter.clone());
  if (not variant_filter) {

    ExecEnv::log().error("OffsetDB::copyFilter; Filter: {} is not a variant_filter.", filter.filterName());

  }

  // Filter the variants.
  std::unique_ptr<OffsetDB> filtered_offset_ptr(std::make_unique<OffsetDB>());
  for (const auto& variant_ptr : variant_vector_) {

    if (variant_filter->applyFilter(*variant_ptr)) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


std::pair<size_t, size_t> kgl::OffsetDB::selfFilter(const BaseFilter &filter) {

  std::shared_ptr<const FilterOffsets> offset_filter = std::dynamic_pointer_cast<const FilterOffsets>(filter.clone());
  if (offset_filter) {

    size_t prior_count = variant_vector_.size();

    auto filtered_offset_ptr = offset_filter->applyFilter(*this);
    variant_vector_ = std::move(filtered_offset_ptr->variant_vector_);

    size_t post_count = variant_vector_.size();

    return { prior_count, post_count};

  }

  std::pair<size_t, size_t> filter_count{0, 0};
  // Only variant filters should be at this level, so we check for this.
  std::shared_ptr<const FilterVariants> variant_filter = std::dynamic_pointer_cast<const FilterVariants>(filter.clone());
  filter_count.first = variant_vector_.size();
  if (not variant_filter) {

    ExecEnv::log().error("OffsetDB::selfFilter; Filter: {} is not a variant_filter.", filter.filterName());
    filter_count.second = filter_count.first;
    return filter_count;

  }
  // Filter the variants.
  OffsetDBArray filtered_variants;
  for (const auto& variant_ptr : variant_vector_) {

    if (variant_filter->applyFilter(*variant_ptr)) {

      filtered_variants.push_back(variant_ptr);

    }

  }
  filter_count.second = filtered_variants.size();
  variant_vector_ = std::move(filtered_variants);

  return filter_count;

}


// setIntersection returns an offset that contains unique variants present in both offsets.
// The VariantEquality flag determines whether variant phase is used in the equality.
std::unique_ptr<kgl::OffsetDB> kgl::OffsetDB::setIntersection( const OffsetDB& offset_B, VariantEquality variant_equality) const {

  std::unique_ptr<OffsetDB> intersection(std::make_unique<OffsetDB>());
  // Generate a set of unique variants from this offset.
  std::unordered_set<std::string> unique_hash_set;
  for (auto const &variant_ptr : getVariantArray()) {

    unique_hash_set.emplace(variant_ptr->equalityHash(variant_equality));

  }

  // Add unique variants from the second offset.
  for (auto const &variant_ptr : offset_B.getVariantArray()) {

    auto find_iter = unique_hash_set.find(variant_ptr->equalityHash(variant_equality));
    if (find_iter != unique_hash_set.end()) {

      intersection->addVariant(variant_ptr);

    }

  }

  return intersection;

}

// setComplement returns a offset that contains variants present in this contig but not present in the complement_contig.
std::unique_ptr<kgl::OffsetDB> kgl::OffsetDB::setComplement( const OffsetDB& offset_B, VariantEquality variant_equality) const {

  std::unique_ptr<OffsetDB> complement(std::make_unique<OffsetDB>());
  // Create a set a complements.
  std::unordered_set<std::string> unique_hash_set;
  for (auto const &variant_ptr : offset_B.getVariantArray()) {

    unique_hash_set.emplace(variant_ptr->equalityHash(variant_equality));

  }

  for (auto const &variant_ptr : getVariantArray()) {

    auto result = unique_hash_set.find(variant_ptr->equalityHash(variant_equality));
    if (result == unique_hash_set.end()) {

      complement->addVariant(variant_ptr);

    }

  }

  return complement;

}
// setUnion returns a offset that contains the set union of unique variants present in this offset and the union_offset.
// Specifically, identical variants, defined by the VariantEquality flag as including phase or not, are ignored.
std::unique_ptr<kgl::OffsetDB> kgl::OffsetDB::setUnion( const OffsetDB& offset_B, VariantEquality variant_equality) const {

  std::unique_ptr<OffsetDB> union_offset(std::make_unique<OffsetDB>());
  // Generate a set of unique variants from this offset.
  std::unordered_set<std::string> unique_hash_set;
  for (auto const &variant_ptr : getVariantArray()) {

    auto[iter, result] = unique_hash_set.emplace(variant_ptr->equalityHash(variant_equality));
    if (result) {

      union_offset->addVariant(variant_ptr);

    }

  }

  // Add unique variants from the second offset.
  for (auto const &variant_ptr : offset_B.getVariantArray()) {

    auto[iter, result] = unique_hash_set.emplace(variant_ptr->equalityHash(variant_equality));
    if (result) {

      union_offset->addVariant(variant_ptr);

    }

  }

  return union_offset;

}

