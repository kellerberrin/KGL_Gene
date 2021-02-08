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


// Ensure max 2 variants per offset.
std::pair<size_t, size_t> kgl::OffsetDB::inSituDiploid() {

  constexpr const static size_t DIPLOID_COUNT{2};

  if (variant_vector_.size() <= DIPLOID_COUNT) {

    return { variant_vector_.size(), variant_vector_.size() };

  }


  std::pair<size_t, size_t> offset_count{0, 0};
  offset_count.first = variant_vector_.size();

  struct DiploidRecord {

    std::vector<std::shared_ptr<const Variant>> variant_vector;
    double frequency{0.0};

  };
  std::unordered_map<std::string, DiploidRecord> hash_map;
  // Construct the hash map.
  for (auto const& variant_ptr : variant_vector_) {

    std::string allele_hash = variant_ptr->alleleHash();
    auto result = hash_map.find(allele_hash);
    if (result != hash_map.end()) {

      auto& [hash, record] = *result;
      record.variant_vector.push_back(variant_ptr);

    } else {

      DiploidRecord diploid_record;
      diploid_record.variant_vector.push_back(variant_ptr);
      auto frequency_opt = FrequencyDatabaseRead::processSuperPopField(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
      if (frequency_opt) {

        diploid_record.frequency = frequency_opt.value();

      } else {

        ExecEnv::log().error("OffsetDB::inSituDiploid; unable to retrieve AF frequency for variant: {}",
                             variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

      }

      auto const [insert_iter, insert_result] = hash_map.try_emplace(allele_hash, diploid_record);
      if (not insert_result) {

        ExecEnv::log().error("OffsetDB::inSituDiploid; could not insert variant hash: {}", allele_hash);

      }

    }

  }

  if (hash_map.empty()) {

    ExecEnv::log().error("OffsetDB::inSituDiploid; unexpected empty hash map");
    offset_count.second = 0;
    return offset_count;

  }

  // Create a map ordered by variant frequency
  std::map<double, DiploidRecord> freq_map;

  for (auto const& [hash, record] : hash_map) {

    freq_map[record.frequency] = record;

  }

  // The most likely highest frequency alleles are last
  auto& [freq, record] = *freq_map.rbegin();
  auto& filtered_variants = record.variant_vector;
  if (filtered_variants.size() > DIPLOID_COUNT) {

    ExecEnv::log().error("OffsetDB::inSituDiploid; Unexpected filtered variant size: {}", record.variant_vector.size());
    for (auto const& variant_ptr : filtered_variants) {

      ExecEnv::log().warn("OffsetDB::inSituDiploid; Unexpected size variant: {}",
                          variant_ptr->output(',', VariantOutputIndex::START_0_BASED, false));

    }

    filtered_variants.resize(DIPLOID_COUNT);

  }

  // Check phase.
  if (filtered_variants.size() == DIPLOID_COUNT) {

    if (filtered_variants.front()->phaseId() == filtered_variants.back()->phaseId()
        and filtered_variants.front()->phaseId() != VariantPhase::UNPHASED) {

      ExecEnv::log().error("OffsetDB::inSituDiploid; Invalid Diploid phase, {}, {}",
                           filtered_variants.front()->output(',', VariantOutputIndex::START_0_BASED, false),
                           filtered_variants.back()->output(',', VariantOutputIndex::START_0_BASED, false));

      std::const_pointer_cast<Variant>(filtered_variants.front())->updatePhaseId(VariantPhase::DIPLOID_PHASE_A);
      std::const_pointer_cast<Variant>(filtered_variants.front())->updatePhaseId(VariantPhase::DIPLOID_PHASE_B);

    }

  }

  variant_vector_ = std::move(filtered_variants);
  offset_count.second = variant_vector_.size();

  return offset_count;

}
