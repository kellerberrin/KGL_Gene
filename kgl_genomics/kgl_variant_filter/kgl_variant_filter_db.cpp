//
// Created by kellerberrin on 23/04/23.
//

#include "kgl_variant_filter_db.h"
#include "kel_interval.h"

#include <ranges>

namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters a population to the listed genomes (if they exist).
// Note that this is a shallow copy of the original population.
// Use selfFilter() or deepCopy() to create a permanent population view.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::PopulationDB> kgl::GenomeListFilter::applyFilter(const PopulationDB& population) const {

  std::unique_ptr<PopulationDB> filtered_population_ptr(std::make_unique<PopulationDB>(population.populationId(), population.dataSource()));

  for (auto const& genome_id : genome_set_) {

    if (population.getMap().contains(genome_id)) {

      auto iter = population.getMap().find(genome_id);
      auto const& [id, genome_ptr] = *iter;

      if (not filtered_population_ptr->addGenome(genome_ptr)) {

        ExecEnv::log().warn("GenomeListFilter::contains; Unexpected duplicate genome: {} in population: {}", genome_id, population.populationId());

      }

    }

  }

  return filtered_population_ptr;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contig Region filter - Uses the half open interval convention [Begin, End).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::ContigDB> kgl::ContigRegionFilter::applyFilter(const ContigDB& contig) const {

  std::unique_ptr<ContigDB> contig_ptr(std::make_unique<ContigDB>(contig.contigId()));

  auto const lower_bound = contig.getMap().lower_bound(start_);
  auto const upper_bound = contig.getMap().lower_bound(end_); //  [start, end)

  auto iter = lower_bound;
  while (iter != contig.getMap().end() and iter != upper_bound) {

    auto const& [offset, offset_ptr] = *iter;

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not contig_ptr->addVariant(variant_ptr)) {

        ExecEnv::log().error("ContigRegionFilter::applyFilter; unable to add variant: {} to contig: {}",
                             variant_ptr->HGVS(), contig_ptr->contigId());

      }

    }

    iter = std::ranges::next(iter, 1, contig.getMap().end());

  }

  return contig_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Contig Modify filter - All variants that modify the region [Begin, End).
// Important - this includes upstream indel deletes that extend into the region.
// The offset these deletes is less than start but the delete variants extends into the region.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::ContigDB> kgl::ContigModifyFilter::applyFilter(const ContigDB& contig) const {

  // Decrement the region start by the heuristic margin.
  ContigOffset_t margin_start = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(start_)-UPSTREAM_DELETE_MARGIN);
  // Return all variants in the [start, end) region.
  auto region_contig_ptr = contig.viewFilter(ContigRegionFilter(start_, end_));

  // Set up the region as an interval.
  const OpenRightInterval specified_region(start_, end_);

  auto const lower_bound = contig.getMap().lower_bound(margin_start);
  auto const upper_bound = contig.getMap().lower_bound(start_); //  [start, end)

  auto iter = lower_bound;
  while (iter != contig.getMap().end() and iter != upper_bound) {

    auto const& [offset, offset_ptr] = *iter;

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      auto [extend_offset, extend_size] = variant_ptr->extentOffset();
      const OpenRightInterval variant_interval(extend_offset, extend_offset + extend_size);

      if (variant_interval.intersects(specified_region)) {

        // Check if a delete variant
        if (variant_ptr->variantType() == VariantType::INDEL_DELETE) {

          if (not region_contig_ptr->addVariant(variant_ptr)) {

            ExecEnv::log().error("ContigModifyFilter::applyFilter; unable to add variant: {} to contig: {}",
                                 variant_ptr->HGVS(), region_contig_ptr->contigId());

          } // If successful add variant.

        } else {

          ExecEnv::log().error("ContigModifyFilter::applyFilter; Unexpected, upstream modifying variant is not an indel DELETE: {}", variant_ptr->HGVS());

        } // If delete variant.

      } // If intersects.

    } // For all variants defined for the offset.

    iter = std::ranges::next(iter, 1, contig.getMap().end());

  } // For all offsets in the margin region.

  return region_contig_ptr;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters all offsets where there are exactly 2 identical variants (disregarding phase).
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::HomozygousFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  if (offset.getVariantArray().size() != 2) {

    return filtered_offset_ptr;

  }
  if (offset.getVariantArray().front()->HGVS() == offset.getVariantArray().back()->HGVS()) {

    filtered_offset_ptr->addVariant(offset.getVariantArray().front());
    filtered_offset_ptr->addVariant(offset.getVariantArray().back());
    return filtered_offset_ptr;

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Ensure max 2 variants per offset.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::DiploidFilter::applyFilter(const OffsetDB& offset) const {

  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  if (offset.getVariantArray().size() <= 2) {

    for (auto const& variant_ptr : offset.getVariantArray()) {

      filtered_offset_ptr->addVariant(variant_ptr);

    }

    return filtered_offset_ptr;

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter unique variants disregarding phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::OffsetDB> kgl::UniqueUnphasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants_;
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    auto variant_hash = variant_ptr->HGVS();
    if (not hashed_variants_.contains(variant_hash)) {

      hashed_variants_.insert(variant_hash);
      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Filter only unique variants including phase.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::OffsetDB> kgl::UniquePhasedFilter::applyFilter(const OffsetDB& offset) const {

  std::unordered_set<std::string> hashed_variants_;
  auto filtered_offset_ptr = std::make_unique<OffsetDB>();
  for (auto const &variant_ptr: offset.getVariantArray()) {

    auto variant_hash = variant_ptr->HGVS_Phase();
    if (not hashed_variants_.contains(variant_hash)) {

      hashed_variants_.insert(variant_hash);
      filtered_offset_ptr->addVariant(variant_ptr);

    }

  }

  return filtered_offset_ptr;

}



