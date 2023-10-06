//
// Created by kellerberrin on 26/07/23.
//

#include "kgl_genome_seq/kgl_seq_variant_filter.h"
#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_variant_filter_coding.h"
#include "kel_interval_unsigned.h"

#include <algorithm>


namespace kgl = kellerberrin::genome;


bool kgl::SequenceVariantFilter::getSortedVariants(const std::shared_ptr<const GenomeDB>& ,
                                                   ContigId_t,
                                                   VariantPhase,
                                                   ContigOffset_t,
                                                   ContigOffset_t,
                                                   OffsetVariantMap &) { return true; }




// Returns a map of unique canonical variants
// Also returns the number of multiple variants found at each offset which are filtered to a single variant.
std::tuple<kgl::OffsetVariantMap, size_t, size_t> kgl::SequenceVariantFilter::getCanonicalVariants(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                                                   const OpenRightUnsigned& variant_interval) {

  OffsetVariantMap offset_variant_map;

  // Filter the variants plus a margin of 200 nucleotides at the beginning of the specified region.
  // To allow for any offset() change when the variants are converted to canonical.
  ContigOffset_t lower = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(variant_interval.lower())-NUCLEOTIDE_CANONICAL_MARGIN);
  // Filter to just variants in the region plus offset margins.
  auto region_contig_ptr = contig_ptr->viewFilter(ContigRegionFilter(lower, variant_interval.upper()));
  // Convert to canonical variants.
  auto canonical_contig_ptr = region_contig_ptr->canonicalContig();
  // Filter to just the variants that will modify the specified region [start, end).
  auto modify_contig_ptr = region_contig_ptr->viewFilter(ContigModifyFilter(variant_interval.lower(), variant_interval.upper()));
  // Get the count of variants modifying the region [start, end).

  // Remove multiple variants, that is minor alleles (SNP) that are not homozygous that share the same offset.
  // Note that if an indel occurs at the same offset as an SNP, this is not a problem, since in canonical form ('1MnI' or '1MnD'),
  // the indel actually modifies the next (+1) offset.
  auto unique_contig_ptr = modify_contig_ptr->viewFilter(HomozygousCodingFilter());
  // Finally filter any variants that are deleted by an upstream delete variant.
  auto no_upstream_delete = unique_contig_ptr->viewFilter(ContigUpstreamFilter());
  // Get the number of unique variants.
  size_t modify_count = (modify_contig_ptr->viewFilter(UniqueUnphasedFilter()))->variantCount();
  // Get the number of upstream deleted variants.
  size_t upstream_deleted = unique_contig_ptr->variantCount() - no_upstream_delete->variantCount();

  // Finally move the unique modifying variants to the offset map.
  for (auto const& [offset, offset_ptr] : no_upstream_delete->getMap()) {

    for (auto const& variant_ptr : offset_ptr->getVariantArray()) {

      // Indels are inserted at the next offset (offset+1) as that is where the insert and delete occurs.
      ContigOffset_t insert_offset{0};
      if (variant_ptr->isSNP()) {

        insert_offset = variant_ptr->offset();

      } else {

        insert_offset = variant_ptr->offset() + 1;

      }

      auto const [insert_iter, result] = offset_variant_map.try_emplace(insert_offset, variant_ptr);
      if (not result) {

        auto const& [find_offset, find_variant] = *offset_variant_map.find(insert_offset);
        ExecEnv::log().error("SequenceVariantFilter::getCanonicalVariants; insert fails at offset: {} insert variant: {}, duplicate map variant: {}"
            , insert_offset, variant_ptr->HGVS(),  find_variant->HGVS());

      }

    } // For all variants at each offset.

  } // for all offsets.

  size_t non_unique_count = modify_count - offset_variant_map.size();

  return { offset_variant_map, non_unique_count, upstream_deleted };

}


void kgl::SequenceVariantFilter::canonicalVariants(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                   const OpenRightUnsigned& sequence_interval) {


  auto [interval_map, non_unique_count, upstream_deleted] = SequenceVariantFilter::getCanonicalVariants(contig_ptr, sequence_interval);
  offset_variant_map_ = std::move(interval_map);
  duplicate_variants_ = non_unique_count;
  downstream_delete_ = upstream_deleted;

}


