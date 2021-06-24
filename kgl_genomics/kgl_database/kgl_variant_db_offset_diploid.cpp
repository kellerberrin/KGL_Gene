//
// Created by kellerberrin on 9/2/21.
//


#include "kgl_variant_db_offset.h"
#include "kgl_variant_db_freq.h"
#include <unordered_set>


namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Fairly complex logic to ensure that diploid populations are
// properly phased and contain a maximum of 2 alleles per location.
// The main (only) problem seems to be that a change of base at a specific location can be due to an SNP
// or an indel. The VCF files often record both these competing allele types.
// The code below looks up the frequency of these competing alleles and retains
// the most likely to occur.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Ensure max 2 variants per offset.
std::pair<size_t, size_t> kgl::OffsetDB::inSituDiploid() {

  constexpr const static size_t DIPLOID_COUNT{2};

  std::pair<size_t, size_t> offset_count{0, 0};
  offset_count.first = variant_vector_.size();

  if (variant_vector_.size() < DIPLOID_COUNT) {

    offset_count.second = variant_vector_.size();
    return offset_count;

  } else if (variant_vector_.size() == DIPLOID_COUNT) {

    // Check phasing if phased.
    if (variant_vector_.front()->phaseId() == variant_vector_.back()->phaseId()
        and variant_vector_.front()->phaseId() != VariantPhase::UNPHASED) {

      // We generally have an Indel and an SNP changing the same reference base.
      // Check to see which one is more common (higher frequency) and retain it, delete the other (VCF flaw).
      auto front_opt = FrequencyDatabaseRead::superPopFrequency(*variant_vector_.front(), FrequencyDatabaseRead::SUPER_POP_ALL_);
      auto back_opt = FrequencyDatabaseRead::superPopFrequency(*variant_vector_.back(), FrequencyDatabaseRead::SUPER_POP_ALL_);

      if (front_opt and back_opt) {

        if (front_opt.value() >= back_opt.value()) {

          // Delete the back element
          variant_vector_.pop_back();

        } else {

          // Delete the front element.
          variant_vector_.erase(variant_vector_.begin());

        }

      } else {
        // Complain and retain the front variant.
        ExecEnv::log().warn("OffsetDB::inSituDiploid; unable to retrieve AF frequency for: {} AND {}",
                            variant_vector_.front()->output(',', VariantOutputIndex::START_0_BASED, true),
                            variant_vector_.back()->output(',', VariantOutputIndex::START_0_BASED, true));
        // Delete the back element
        variant_vector_.pop_back();

      }

    }

    offset_count.second = variant_vector_.size();
    return offset_count;

  }

  // More than 2 variants at the location.
  struct DiploidRecord {

    OffsetDBArray variant_vector;
    double frequency{0.0};

  };
  std::unordered_map<std::string, DiploidRecord> hash_map;
  // Construct the hash map.
  for (auto const& variant_ptr : variant_vector_) {

    std::string variant_hash = variant_ptr->variantHash();
    auto result = hash_map.find(variant_hash);
    if (result != hash_map.end()) {

      auto& [hash, record] = *result;
      record.variant_vector.push_back(variant_ptr);

    } else {

      DiploidRecord diploid_record;
      diploid_record.variant_vector.push_back(variant_ptr);
      auto frequency_opt = FrequencyDatabaseRead::superPopFrequency(*variant_ptr, FrequencyDatabaseRead::SUPER_POP_ALL_);
      if (frequency_opt) {

        diploid_record.frequency = frequency_opt.value();

      } else {

        ExecEnv::log().warn("OffsetDB::inSituDiploid; unable to retrieve AF frequency for variant: {}",
                            variant_ptr->output(',', VariantOutputIndex::START_0_BASED, true));
        diploid_record.frequency = 0.0;

      }

      auto const [insert_iter, insert_result] = hash_map.try_emplace(variant_hash, diploid_record);
      if (not insert_result) {

        ExecEnv::log().error("OffsetDB::inSituDiploid; could not insert variant hash: {}", variant_hash);

      }

    }

  }


  // Create a map ordered by variant frequency
  std::map<double, DiploidRecord> freq_map;
  for (auto const& [hash, record] : hash_map) {

    freq_map[record.frequency] = record;

  }

  if (freq_map.empty()) {

    ExecEnv::log().warn("OffsetDB::inSituDiploid; unexpected empty frequency map");
    variant_vector_.clear();
    offset_count.second = variant_vector_.size();
    return offset_count;

  }

  // The most likely highest frequency alleles are last
  auto const& [freq, record] = *freq_map.rbegin();
  OffsetDBArray filtered_variants = record.variant_vector;

  // If unphased just resize
  if (filtered_variants.front()->phaseId() == VariantPhase::UNPHASED
      and filtered_variants.size() > DIPLOID_COUNT) {

    filtered_variants.resize(DIPLOID_COUNT);

  } else if (filtered_variants.front()->phaseId() != VariantPhase::UNPHASED) {

    // look for complementary phases.
    bool first_pass = true;
    std::string first_hash = filtered_variants.front()->variantPhaseHash();
    OffsetDBArray phase_subset;
    for (auto& variant_ptr :  filtered_variants) {

      if (first_pass) {

        phase_subset.push_back(variant_ptr);
        first_pass = false;
        continue;

      }

      if (variant_ptr->variantPhaseHash() != first_hash) {

        phase_subset.push_back(variant_ptr);
        break;

      }

    }

    filtered_variants = phase_subset;

  }

  // If 2 variants, check phase
  if (filtered_variants.size() == DIPLOID_COUNT) {

    if (filtered_variants.front()->phaseId() == filtered_variants.back()->phaseId()
        and filtered_variants.front()->phaseId() != VariantPhase::UNPHASED) {

      ExecEnv::log().warn("OffsetDB::inSituDiploid; Invalid Diploid phase, {}, {}",
                          filtered_variants.front()->output(',', VariantOutputIndex::START_0_BASED, true),
                          filtered_variants.back()->output(',', VariantOutputIndex::START_0_BASED, true));

      filtered_variants.clear();

    }

  }

  variant_vector_ = filtered_variants;
  offset_count.second = variant_vector_.size();

  return offset_count;

}
