//
// Created by kellerberrin on 4/2/21.
//


#include "kgl_variant_db_offset.h"
#include "kgl_variant_db_freq.h"
#include <unordered_set>


namespace kgl = kellerberrin::genome;



std::pair<size_t, size_t> kgl::OffsetDB::inSituFilter(const VariantFilter &filter) {

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

    case FilterType::DiploidFilter:
      return inSituDiploid();

    case FilterType::UniqueUnphasedFilter:
    case FilterType::UniquePhasedFilter:
      return inSituUnique(filter);

    case FilterType::TrueFilter:
      return  {variant_vector_.size(), variant_vector_.size()};

    case FilterType::FalseFilter:
      variant_vector_.clear();
      return {variant_vector_.size(), 0};

    case FilterType::GenomeFilter:   // Not implemented at this level.
      return {variant_vector_.size(), variant_vector_.size()};

    default:
      ExecEnv::log().error("OffsetDB::inSituFilter; Unknown filter type");
      return {variant_vector_.size(), variant_vector_.size()};

  }

}

std::pair<size_t, size_t> kgl::OffsetDB::inSituGeneral(const VariantFilter &filter) {

  std::pair<size_t, size_t> offset_count{0, 0};

  offset_count.first = variant_vector_.size();
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


std::pair<size_t, size_t> kgl::OffsetDB::inSituUnique(const VariantFilter &filter) {

  std::pair<size_t, size_t> offset_count{0, 0};
  std::unordered_set<std::string> hash_set;

  offset_count.first = variant_vector_.size();
  OffsetDBArray filtered_variants;
  filtered_variants.reserve(offset_count.first);
  bool unique_phased = filter.filterType() == FilterType::UniquePhasedFilter;

  for (auto const &variant_ptr : variant_vector_) {

    auto variant_hash = unique_phased ? variant_ptr->variantPhaseHash() : variant_ptr->variantHash();
    auto const &result = hash_set.find(variant_hash);
    if (result == hash_set.end()) {

      filtered_variants.push_back(variant_ptr);

      auto const&[iter, insert_result] = hash_set.emplace(std::move(variant_hash));
      if (not insert_result) {

        ExecEnv::log().error("OffsetDB::inSituUnique; variant hash: {} cannot be inserted",
                             (unique_phased ? variant_ptr->variantPhaseHash() : variant_ptr->variantHash()));

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

