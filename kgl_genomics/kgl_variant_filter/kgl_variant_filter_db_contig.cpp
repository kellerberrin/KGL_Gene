//
// Created by kellerberrin on 11/08/23.
//

#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kel_interval.h"


namespace kgl = kellerberrin::genome;



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


  // Return all variants in the [start, end) region.
  auto region_contig_ptr = contig.viewFilter(ContigRegionFilter(start_, end_));

  // Set up the region as an interval.
  const OpenRightInterval specified_region(start_, end_);

  // Decrement the region start by the heuristic margin.
  ContigOffset_t margin_start = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(start_)-UPSTREAM_DELETE_MARGIN);
  auto const lower_bound = contig.getMap().lower_bound(margin_start);

  ContigOffset_t margin_end = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(start_)-1);
  auto const upper_bound = contig.getMap().lower_bound(margin_end); //  Search [start-margin, start-1)

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

          ExecEnv::log().error("ContigModifyFilter::applyFilter; Unexpected, upstream of [ {}, {}) modifying variant is not an indel DELETE: {}"
                               , start_, end_, variant_ptr->HGVS());

        } // If delete variant.

      } // If intersects.

    } // For all variants defined for the offset.

    iter = std::ranges::next(iter, 1, contig.getMap().end());

  } // For all offsets in the margin region.

  return region_contig_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns a contig containing all variants in *this contig that match the template contig.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This search algorithm has n^2 complexity.
// The variants in the template contig are unique. Variant phase is disregarded.
std::unique_ptr<kgl::ContigDB> kgl::ContigTemplateFilter::applyFilter(const ContigDB& contig) const {

  auto found_contig_ptr = std::make_unique<ContigDB>(reference_ptr_->contigId());

  for (auto const& [offset, offset_ptr] : reference_ptr_->getMap()) {

    auto find_iter = contig.getMap().find(offset);
    if (find_iter != contig.getMap().end()) {

      // Create a set of allele hashs to search.
      std::unordered_set<std::string> search_hash;
      for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

        search_hash.insert(variant_ptr->HGVS());

      }

      // Search the set of hashs.
      auto const& [this_offset, this_offset_ptr] = *find_iter;
      for (auto const& this_variant_ptr : this_offset_ptr->getVariantArray()) {

        if (search_hash.contains(this_variant_ptr->HGVS())) {

          if (not found_contig_ptr->addVariant(this_variant_ptr)) {

            ExecEnv::log().error( "ContigDB::findContig; cannot add variant: {}", this_variant_ptr->HGVS());

          }

        }

      }

    } // if this offset

  } // for all template offset


  return found_contig_ptr;

}

