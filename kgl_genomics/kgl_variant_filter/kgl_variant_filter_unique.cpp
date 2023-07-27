//
// Created by kellerberrin on 17/07/23.
//

#include "kgl_variant_filter_unique.h"
#include "kgl_variant_filter_info.h"

#include <ranges>


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Selection logic is somewhat involved because canonical INDEL variants actually operate on the NEXT (offset+1) offset.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::ContigDB> kgl::RandomUniqueFilter::filterUnique(const ContigDB &contig) const {

  std::unique_ptr<ContigDB> filtered_contig_ptr = std::make_unique<ContigDB>(contig.contigId());

  // Examine all offsets
  std::vector<std::shared_ptr<const Variant>> current_offset_vector;
  std::vector<std::shared_ptr<const Variant>> next_offset_vector;

  auto contig_iter = contig.getMap().begin();
  while (contig_iter != contig.getMap().end()) {

    current_offset_vector = next_offset_vector;
    next_offset_vector.clear();

    auto const& [current_offset, current_offset_ptr] = *contig_iter;

    bool next_offset_active{false};
    auto next_iter = std::ranges::next(contig_iter, 1, contig.getMap().end());
    if (next_iter != contig.getMap().end()) {

      auto const& [next_offset, next_offset_ptr] = *next_iter;
      next_offset_active = (current_offset + 1) == next_offset;

    }

    for (auto const& variant_ptr : current_offset_ptr->getVariantArray()) {

      if (not variant_ptr->isCanonical()) {

        ExecEnv::log().error("UniqueOffsetFilter::applyFilter; variant NOT canonical: {}", variant_ptr->HGVS());
        continue;

      }

      if (variant_ptr->isSNP()) {

        current_offset_vector.push_back(variant_ptr);

      } else if (next_offset_active) {

        next_offset_vector.push_back(variant_ptr);

      }

    } // for offset variants.

    // Select the candidate variants by allele frequency.
    if (not current_offset_vector.empty()) {

      auto selected_variant = selectUnique(current_offset_vector);

      if (not filtered_contig_ptr->addVariant(selected_variant)) {

        ExecEnv::log().error("UniqueOffsetFilter::applyFilter; unable to add variant: {} to contig: {}",
                             selected_variant->HGVS(), filtered_contig_ptr->contigId());

      }

    }

    // Increment to next offset.
    contig_iter = std::ranges::next(contig_iter, 1, contig.getMap().end());

  }

  return filtered_contig_ptr;

}


std::shared_ptr<const kgl::Variant> kgl::RandomUniqueFilter::selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("UniqueOffsetFilter::selectRandom; selection vector is empty - cannot continue");

  }

  auto selected_variant = variant_vector.front();

  return selected_variant;

}


std::shared_ptr<const kgl::Variant> kgl::FrequencyUniqueFilter::selectUnique(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("UniqueOffsetFilter::selectRandom; selection vector is empty - cannot continue");

  }

  std::multimap<double, std::shared_ptr<const kgl::Variant>> frequency_map;
  for (auto const& variant_ptr : variant_vector) {

    double frequency = getFrequency(variant_ptr);
    frequency_map.insert({frequency, variant_ptr});

  }

  auto [selected_frequency, selected_variant] = *frequency_map.rbegin();

  return selected_variant;

}


double kgl::FrequencyUniqueFilter::getFrequency(const std::shared_ptr<const Variant>& variant_ptr) const {


  size_t alt_count = variant_ptr->evidence().altVariantCount();
  size_t alt_index = variant_ptr->evidence().altVariantIndex();

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<std::vector<double>>(*variant_ptr, AF_FIELD_);
  if (info_opt) {

    std::vector<double> info_vector = std::move(info_opt.value());
    if (info_vector.size() != alt_count) {

      ExecEnv::log().error("FrequencyUniqueFilter::getFrequency; AF vector size: {}, not equal alt variant count: {}, Info field: {}",
                           info_vector.size(), alt_count, AF_FIELD_);

      return 0.0;

    }
    if (info_vector.size() <= alt_index) {

      ExecEnv::log().error("FrequencyUniqueFilter::getFrequency; alt variant index: {} out of range for vector size:{}, Info field: {}",
                           alt_index, info_vector.size(), AF_FIELD_);

    }

    return info_vector[alt_index];

  }

  return 0.0;

}
