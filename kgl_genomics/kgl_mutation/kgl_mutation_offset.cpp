//
// Created by kellerberrin on 26/07/23.
//

#include "kgl_mutation_offset.h"
#include "kgl_variant_filter_db.h"
#include "kgl_variant_filter_unique.h"
#include "kel_interval.h"

#include <algorithm>


namespace kgl = kellerberrin::genome;


bool kgl::MutationOffset::getSortedVariants( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                             ContigId_t contig_id,
                                             VariantPhase,
                                             ContigOffset_t start,
                                             ContigOffset_t end,
                                             OffsetVariantMap &variant_map) {


  auto find_iter = genome_ptr->getMap().find(contig_id);

  if (find_iter == genome_ptr->getMap().end()) {

    ExecEnv::log().error("MutationOffset::getSortedVariants; Contig Id: {} not found in Genome Variant: {}", contig_id, genome_ptr->genomeId());
    return false;

  }

  std::shared_ptr<ContigDB> contig_ptr = find_iter->second;

  auto [contig_map, result] = MutationOffset::getCanonicalVariants(contig_ptr, start, end);

  variant_map = std::move(contig_map);

  return result;

}


std::pair<kgl::OffsetVariantMap, bool> kgl::MutationOffset::getCanonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                                  ContigOffset_t start,
                                                                                  ContigOffset_t end) {

  OffsetVariantMap offset_variant_map;

  // Filter the variants plus a margin of 200 nucleotides at the beginning of the specified region.
  // To allow for any offset() change when the variants are converted to canonical.
  size_t lower = std::max<int64_t>(0, static_cast<int64_t>(start)-NUCLEOTIDE_CANONICAL_MARGIN);
  // Filter to just variants in the region plus offset margins.
  auto region_contig_ptr = contig_ptr->viewFilter(ContigRegionFilter(lower, end));
  // Remove any duplicate (homozygous) identical variants.
  auto no_homo_contig_ptr = region_contig_ptr->viewFilter(UniqueUnphasedFilter());
  // Convert to canonical variants.
  auto canonical_contig_ptr = no_homo_contig_ptr->canonicalContig();
  // Remove multiple variants, that is minor alleles (SNP) that are not homozygous that share the same offset.
  // Note that if an indel occurs at the same offset as an SNP, this is not a problem, since in canonical form ('1MnI' or '1MnD'),
  // the indel actually modifies the next (+1) offset.
  std::shared_ptr<const ContigDB> unique_contig_ptr = canonical_contig_ptr->viewFilter(FrequencyUniqueFilter());
  // Finally remove any excess variants that do not intersect with [start, end). This is simple for SNP and inserts.
  // However deletes can extend from before the start of the interval.
  OpenRightInterval variant_region(start, end);
  for (auto const& [offset, offset_ptr] : unique_contig_ptr->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      auto const [offset_extent, extent_size] = variant_ptr->extentOffset();
      OpenRightInterval variant_footprint(offset_extent, offset_extent + extent_size);
      // If the variant modifies the region [start,end)
      if (variant_region.intersects(variant_footprint)) {

        // Indels are inserted at the next offset (offset+1) as that is where the insert and delete occurs.
        ContigOffset_t insert_offset{0};
        if (variant_ptr->isSNP()) {

          insert_offset = variant_ptr->offset();

        } else {

          insert_offset = variant_ptr->offset() + 1;

        }

        auto const [insert_iter, result] = offset_variant_map.try_emplace(insert_offset, variant_ptr);
        if (not result) {

          ExecEnv::log().error("MutationOffset::getCanonicalVariants; duplicate variant found: {}", variant_ptr->HGVS());

        }

      } // Variant intersects [start, end)

    } // For all variants at each offset.

  } // for all offsets.

  return { offset_variant_map, true };

}

