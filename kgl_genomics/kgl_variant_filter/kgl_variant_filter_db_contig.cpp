//
// Created by kellerberrin on 11/08/23.
//

#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kel_interval_unsigned.h"

#include <ranges>


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Contig Region filter - Uses the half open interval convention [Begin, End).
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<kgl::ContigDB> kgl::ContigRegionFilter::applyFilter(const ContigDB& contig) const {

  std::unique_ptr<ContigDB> contig_ptr(std::make_unique<ContigDB>(contig.contigId()));

  auto const lower_bound = contig.getMap().lower_bound(start_);
  auto const upper_bound = contig.getMap().lower_bound(end_); //  [start, end)

  for (auto const& [offset, offset_ptr] : std::ranges::subrange(lower_bound, upper_bound)) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      if (not contig_ptr->addVariant(variant_ptr)) {

        ExecEnv::log().error("ContigRegionFilter::applyFilter; unable to add variant: {} to contig_ref_ptr: {}",
                             variant_ptr->HGVS(), contig_ptr->contigId());

      }

    }

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

  std::unique_ptr<ContigDB> region_contig_ptr(std::make_unique<ContigDB>(contig.contigId()));

  OpenRightUnsigned specified_region{start_, end_};

  // Decrement the region start by the heuristic margin.
  ContigOffset_t margin_start = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(start_)-UPSTREAM_DOWNSTREAM_MARGIN_);
  ContigOffset_t margin_end = end_ + UPSTREAM_DOWNSTREAM_MARGIN_;

  auto const lower_bound = contig.getMap().lower_bound(margin_start);
  auto const upper_bound = contig.getMap().lower_bound(margin_end); //  Search [start-margin, start-1)

  for (auto const& [offset, offset_ptr] : std::ranges::subrange(lower_bound, upper_bound)) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      // Note that memberInterval() is not the same as modifyInterval().
      // This function is used to determine if a variant actually modifies the interval of interest
      // rather than just modifying a region adjacent (Insert) to the interval and translating it's offset.
      auto [variant_type, variant_interval] = variant_ptr->memberInterval();

      bool add_variant{false};
      switch (variant_type) {

        // Delete indels only need to intersect the region of interest to modify it.
        case VariantType::INDEL_DELETE:
          add_variant = specified_region.intersects(variant_interval);
          break;

        case VariantType::SNP:
        case VariantType::INDEL_INSERT:
          add_variant = specified_region.containsInterval(variant_interval);
          break;

      }

      if (add_variant) {

        if (not region_contig_ptr->addVariant(variant_ptr)) {

          ExecEnv::log().error("ContigModifyFilter::applyFilter; unable to add variant: {} to contig_ref_ptr: {}",
                               variant_ptr->HGVS(), region_contig_ptr->contigId());

        } // If successful add variant.

      } // If add_variant flag is set.

    } // For all variants defined for the offset.

  } // For all offsets in the margin region.

  return region_contig_ptr;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter out all variants that will be deleted by an upstream Delete variant.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::unique_ptr<kgl::ContigDB> kgl::ContigUpstreamFilter::applyFilter(const ContigDB& contig) const {

  std::unique_ptr<ContigDB> contig_ptr(std::make_unique<ContigDB>(contig.contigId()));

  for (auto const& [offset, offset_ptr] : contig.getMap()) {

    for (auto const &variant_ptr: offset_ptr->getVariantArray()) {

      auto const [variant_type, member_interval] = variant_ptr->memberInterval();

      OpenRightUnsigned lower_bound_key{member_interval.lower(), member_interval.lower()}; // Zero sized.
      auto const lower_bound = upstream_delete_map_.lower_bound(lower_bound_key);
      auto const upper_bound = upstream_delete_map_.end();

      bool upstream_delete{false};

      for (auto const& [delete_interval, delete_variant_ptr] : std::ranges::subrange(lower_bound, upper_bound)) {

        if (delete_interval.intersects(member_interval)) {

          upstream_delete = true;
          break;

        }

      }

      if (not upstream_delete) {

        if (not contig_ptr->addVariant(variant_ptr)) {

          ExecEnv::log().error("ContigUpstreamFilter::applyFilter; unable to add variant: {} to contig_ref_ptr: {}",
                               variant_ptr->HGVS(), contig_ptr->contigId());

        }

        if (variant_type == VariantType::INDEL_DELETE) {

          upstream_delete_map_.emplace(member_interval, variant_ptr);

        } // If delete.

      } // If not upstream delete.

    } // For all variants in the offset.

  } // For all Offsets.

  return contig_ptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns a contig containing all variants in *this contig that match the template contig_ref_ptr.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This search algorithm has n^2 complexity.
// The variants in the template contig_ref_ptr are unique. Variant phase is disregarded.
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

