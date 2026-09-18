//
// Created by kellerberrin on 17/07/23.
//

#include "kgl_variant_filter_coding.h"
#include "kgl_variant_filter_info.h"
#include "kgl_variant_filter_db_offset.h"


#include <ranges>
#include <utility>


namespace kgl = kellerberrin::genome;


/// Filter unique variants for each offset.
///
/// Selection logic is somewhat convoluted because canonical INDEL variants actually operate on the NEXT (offset+1) offset.
/// Indels at offset N are selected together with the variants at N+1; pending indels are flushed on offset gaps > 1
/// and at loop end.
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
        current_offset_vector = std::move(indel_offset_vector);

      } else {
        // Else select unique
        contigVector(filtered_contig_ptr, indel_offset_vector);

      }

      indel_offset_vector.clear();

    }

    for (auto const& variant_ptr : current_offset_ptr->getVariantArray()) {

      if (not variant_ptr->isCanonical()) {

        ExecEnv::log().error("RandomCodingFilter::filterUnique; variant NOT canonical: {}", variant_ptr->HGVS());
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

    // A single candidate is taken as is; multiple candidates are resolved by the selectUnique() strategy.
    auto selected_variant = offset_vector.size() == 1 ? offset_vector.front() : selectUnique(offset_vector);

    // Add to the filtered contig object.
    if (not filtered_contig_ptr->addVariant(selected_variant)) {

      ExecEnv::log().error("RandomCodingFilter::contigVector; unable to add variant: {} to contig_ref_ptr: {}",
                           selected_variant->HGVS(), filtered_contig_ptr->contigId());

    }

  }

}


/// Deterministically selects the first candidate; the class name is historical.
std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectRandom(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("RandomCodingFilter::selectRandom; selection vector is empty - cannot continue");

  }

  // Deterministically selects the first candidate (the 'random' class name is historical).
  return variant_vector.front();

}


/// Returns the allele frequency of the variant's alt allele from the AF info field (0.0 if absent or malformed).
double kgl::RandomCodingFilter::getFrequency(const std::shared_ptr<const Variant>& variant_ptr) {

  size_t alt_count = variant_ptr->evidence().altVariantCount();
  size_t alt_index = variant_ptr->evidence().altVariantIndex();

  auto info_opt = InfoEvidenceAnalysis::getTypedInfoData<std::vector<double>>(*variant_ptr, AF_FIELD_);
  if (info_opt) {

    std::vector<double> info_vector = std::move(info_opt.value());
    if (info_vector.size() != alt_count) {

      ExecEnv::log().error("RandomCodingFilter::getFrequency; AF vector size: {}, not equal alt variant count: {}, Info field: {}",
                           info_vector.size(), alt_count, AF_FIELD_);

      return 0.0;

    }
    if (info_vector.size() <= alt_index) {

      ExecEnv::log().error("RandomCodingFilter::getFrequency; alt variant index: {} out of range for vector size:{}, Info field: {}",
                           alt_index, info_vector.size(), AF_FIELD_);

      // Previously fell through to an out-of-bounds vector access; return 0.0 (lowest frequency) instead.
      return 0.0;

    }

    return info_vector[alt_index];

  }

  return 0.0;

}

/// Selects the most frequently occurring variant (highest probability of occurring).
std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectFrequency(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("RandomCodingFilter::selectFrequency; selection vector is empty - cannot continue");

  }

  // A multimap is used deliberately: rbegin() selects the LAST inserted among equal frequencies,
  // which is the historical tie-break behaviour (it can select a different phase of an otherwise
  // identical variant, which feeds PhaseFilter/HeterozygousFilter downstream).
  std::multimap<double, std::shared_ptr<const kgl::Variant>> frequency_map;
  for (auto const& variant_ptr : variant_vector) {

    double frequency = getFrequency(variant_ptr);
    frequency_map.insert({frequency, variant_ptr});

  }

  auto const& [selected_frequency, selected_variant] = *frequency_map.rbegin();

  return selected_variant;

}


/// Selects among deduplicated candidates by frequency.
/// Note: homozygous preference is currently NOT implemented; the selector is frequency-only among
/// the phase-deduplicated candidates (a homozygous pair collapses to a single variant).
std::shared_ptr<const kgl::Variant> kgl::RandomCodingFilter::selectHomozygous(const std::vector<std::shared_ptr<const Variant>>& variant_vector) const {

  if (variant_vector.empty()) {

    ExecEnv::log().critical("RandomCodingFilter::selectHomozygous; selection vector is empty - cannot continue");

  }

  OffsetDB offset_db;
  for (auto const& variant_ptr : variant_vector) {

    offset_db.addVariant(variant_ptr);

  }

  auto unique_variants = offset_db.viewFilter(UniqueUnphasedFilter());

  if (not unique_variants->getVariantArray().empty()) {

    // getVariantArray() returns the exact const vector<shared_ptr<const Variant>>& that selectFrequency accepts.
    return selectFrequency(unique_variants->getVariantArray());

  }

  return selectFrequency(variant_vector);

}
