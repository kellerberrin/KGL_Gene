//
// Created by kellerberrin on 4/2/21.
//


#include "kgl_variant_db_offset.h"
#include "kgl_variant_db_freq.h"
#include <unordered_set>


namespace kgl = kellerberrin::genome;



std::pair<size_t, size_t> kgl::OffsetDB::selfFilter(const BaseFilter &filter) {

  switch(filter.filterType()) {

    case FilterType::InfoFilter:
    case FilterType::VepSubStringFilter:
    case FilterType::InfoGEQIntegerFilter:
    case FilterType::InfoGEQFloatFilter:
    case FilterType::InfoSubStringFilter:
    case FilterType::InfoBooleanFilter:
    case FilterType::RefAltCountFilter:
    case FilterType::DPCountFilter:
    case FilterType::PhaseFilter:
    case FilterType::PassFilter:
    case FilterType::SNPFilter:
    case FilterType::ContigFilter:
    case FilterType::RegionFilter:
    case FilterType::NotFilter:
    case FilterType::AndFilter:
    case FilterType::OrFilter:
      return inSituGeneral(filter);

    case FilterType::HomozygousFilter:
      return inSituHomozygous();

    case FilterType::DiploidFilter:
      return inSituDiploid();

    case FilterType::UniqueUnphasedFilter:
    case FilterType::UniquePhasedFilter:
      return inSituUnique(filter);

    case FilterType::TrueFilter:
      return  {variant_vector_.size(), variant_vector_.size()};

    case FilterType::FalseFilter: {

      size_t vector_size = variant_vector_.size();
      variant_vector_.clear();
      return { vector_size, 0};

    }

    case FilterType::GenomeFilter:   // Not implemented at this level.
      return {variant_vector_.size(), variant_vector_.size()};

    default:
      ExecEnv::log().error("OffsetDB::selfFilter; Unknown filter type");
      return {variant_vector_.size(), variant_vector_.size()};

  }

}

std::shared_ptr<kgl::OffsetDB> kgl::OffsetDB::copyFilter(const BaseFilter &filter) const {

  std::shared_ptr<OffsetDB> filtered_offset_ptr(std::make_shared<OffsetDB>());

  std::shared_ptr<const FilterOffsets> offset_filter = std::dynamic_pointer_cast<const FilterOffsets>(filter.clone());
  if (offset_filter) {

    if (offset_filter->applyFilter(*this)) {

      // Called recursively, copy all variants.
      return copyFilter(TrueFilter());

    }

    // Else return an empty offset object.
    return filtered_offset_ptr;

  }

  // Only variant filters should be at this level, so we check for this.
  std::shared_ptr<const FilterVariants> variant_filter = std::dynamic_pointer_cast<const FilterVariants>(filter.clone());
  if (not variant_filter) {

    ExecEnv::log().error("OffsetDB::copyFilter; Filter: {} is not a variant_filter.", filter.filterName());

  }
  // Filter the variants.
  for (const auto& variant_ptr : variant_vector_) {

    if (variant_filter->applyFilter(*variant_ptr)) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


std::pair<size_t, size_t> kgl::OffsetDB::inSituGeneral(const BaseFilter &filter) {

  std::pair<size_t, size_t> offset_count{0, 0};

  offset_count.first = variant_vector_.size();
  OffsetDBArray filtered_variants;
  filtered_variants.reserve(offset_count.first);

  auto variant_filter_ptr = std::dynamic_pointer_cast<const FilterVariants>(filter.clone());

  for (auto const &variant_ptr : variant_vector_) {

    if (variant_filter_ptr->applyFilter(*variant_ptr)) {

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

// If there are 2 identical variants, disregarding phase, then the variants are retained else they are deleted.
std::pair<size_t, size_t> kgl::OffsetDB::inSituHomozygous() {

  std::pair<size_t, size_t> offset_count{0, 0};

  offset_count.first = variant_vector_.size();

  if (offset_count.first == 2) {

    if (variant_vector_[0]->HGVS() == variant_vector_[0]->HGVS()) {

      offset_count.second = offset_count.first;
      return offset_count;

    }

  }

  variant_vector_.clear();
  offset_count.second = 0;

  return offset_count;

}


std::pair<size_t, size_t> kgl::OffsetDB::inSituUnique(const BaseFilter &filter) {

  std::pair<size_t, size_t> offset_count{0, 0};
  std::unordered_set<std::string> hash_set;

  offset_count.first = variant_vector_.size();
  OffsetDBArray filtered_variants;
  filtered_variants.reserve(offset_count.first);
  bool unique_phased = filter.filterType() == FilterType::UniquePhasedFilter;

  for (auto const &variant_ptr : variant_vector_) {

    auto variant_hash = unique_phased ? variant_ptr->HGVS_Phase() : variant_ptr->HGVS();
    auto const &result = hash_set.find(variant_hash);
    if (result == hash_set.end()) {

      filtered_variants.push_back(variant_ptr);

      auto const&[iter, insert_result] = hash_set.emplace(std::move(variant_hash));
      if (not insert_result) {

        ExecEnv::log().error("OffsetDB::inSituUnique; variant hash: {} cannot be inserted",
                             (unique_phased ? variant_ptr->HGVS_Phase() : variant_ptr->HGVS()));

      }

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

