//
// Created by kellerberrin on 17/07/23.
//

#include "kgl_variant_filter_coding.h"
#include "kgl_variant_filter_info.h"
#include "kgl_variant_filter_db_offset.h"


#include <ranges>


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Selection logic is somewhat convoluted because canonical INDEL variants actually operate on the NEXT (offset+1) offset.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::ContigDB> kgl::RandomCodingFilter::filterUnique(const ContigDB &contig) const {

  std::unique_ptr<ContigDB> filtered_contig_ptr = std::make_unique<ContigDB>(contig.contigId());

  // Examine all offsets
  std::vector<std::shared_ptr<const Variant>> current_offset_vector;
  std::vector<std::shared_ptr<const Variant>> indel_offset_vector;

  for (auto const& [current_offset, current_offset_ptr] : contig.getMap()) {

    current_offset_vector.clear();
    if (not indel_offset_vector.empty()) {

      if (current_offset == indel_offset_vector.front()->offset() + 1) {
        // Next +1 offset add to current
        current_offset_vector = indel_offset_vector;

      } else {
        // Else select unique
        contigVector(filtered_contig_ptr, indel_offset_vector);

      }

      indel_offset_vector.clear();

    }

    for (auto const& variant_ptr : current_offset_ptr->getVariantArray()) {

      if (not variant_ptr->isCanonical()) {

        ExecEnv::log().error("UniqueOffsetFilter::applyFilter; variant NOT canonical: {}", variant_ptr->HGVS());
        continue;

      }

      if (variant_ptr->isSNP()) {

        current_offset_vector.push_back(variant_ptr);

      } else { //indel

        indel_offset_vector.push_back(variant_ptr);

      }

    } // For all offset variants.

    contigVector(filtered_contig_ptr, current_offset_vector);

  } // For all offsets.

  contigVector(filtered_contig_ptr, indel_offset_vector);

  return filtered_contig_ptr;

}

void kgl::RandomCodingFilter::contigVector(std::unique_ptr<ContigDB>& filtered_contig_ptr,
                                           std::vector<std::shared_ptr<const Variant>>& offset_vector) const {


  // Select the candidate variants by allele frequency.
  if (not offset_vector.empty()) {

    if (offset_vector.size() == 1) {

      // Don't need to select a single variant.
      auto selected_variant = offset_vector.front();
      // Add to the filtered contig_ref_ptr object.
      if (not filtered_contig_ptr->addVariant(selected_variant)) {

        ExecEnv::log().error("UniqueOffsetFilter::applyFilter; unable to add variant: {} to contig_ref_ptr: {}",
                             selected_variant->HGVS(), filtered_contig_ptr->contigId());

      }

    } else {

      // Else select by frequency.
      auto selected_variant = selectUnique(offset_vector);
      // Add to the filtered contig_ref_ptr object.
      if (not filtered_contig_ptr->addVariant(selected_variant)) {

        ExecEnv::log().error("UniqueOffsetFilter::applyFilter; unable to add variant: {} to contig_ref_ptr: {}",
                             selected_variant->HGVS(), filtered_contig_ptr->contigId());

      }

    }

  }

}


std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectRandom(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("UniqueOffsetFilter::selectRandom; selection vector is empty - cannot continue");

  }

  auto selected_variant = variant_vector.front();

  return selected_variant;

}



double kgl::RandomCodingFilter::getFrequency(const std::shared_ptr<const Variant>& variant_ptr) {


  size_t alt_count = variant_ptr->evidence().altVariantCount();
  size_t alt_index = variant_ptr->evidence().altVariantIndex();

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<std::vector<double>>(*variant_ptr, AF_FIELD_);
  if (info_opt) {

    std::vector<double> info_vector = std::move(info_opt.value());
    if (info_vector.size() != alt_count) {

      ExecEnv::log().error("FrequencyCodingFilter::getFrequency; AF vector size: {}, not equal alt variant count: {}, Info field: {}",
                           info_vector.size(), alt_count, AF_FIELD_);

      return 0.0;

    }
    if (info_vector.size() <= alt_index) {

      ExecEnv::log().error("FrequencyCodingFilter::getFrequency; alt variant index: {} out of range for vector size:{}, Info field: {}",
                           alt_index, info_vector.size(), AF_FIELD_);

    }

    return info_vector[alt_index];

  }

  return 0.0;

}

std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectFrequency(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

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


std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectHomozygous(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("UniqueOffsetFilter::selectRandom; selection vector is empty - cannot continue");

  }

  OffsetDB offset_db;
  for (auto const& variant_ptr : variant_vector) {

    offset_db.addVariant(variant_ptr);

  }

  auto homozygous_offset_ptr = offset_db.viewFilter(HomozygousFilter());
  auto unique_homozygous  = offset_db.viewFilter(UniqueUnphasedFilter());

  if (not unique_homozygous->getVariantArray().empty()) {

    std::vector<std::shared_ptr<const Variant>> homozygous_variant_vector;
    for (auto const& variant_ptr : unique_homozygous->getVariantArray()) {

      homozygous_variant_vector.push_back(variant_ptr);

    }

    return selectFrequency(homozygous_variant_vector);

  }

  return selectFrequency(variant_vector);

}


