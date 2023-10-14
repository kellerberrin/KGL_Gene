//
// Created by kellerberrin on 26/07/23.
//

#include "kgl_mutation_variant_filter.h"
#include "kgl_variant_filter_db_contig.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_variant_filter_db_variant.h"

#include "kgl_variant_filter_coding.h"

#include <algorithm>


namespace kgl = kellerberrin::genome;


// Returns a map of unique canonical variants
// Also returns the number of multiple variants found at each offset which are filtered to a single variant.
std::pair<kgl::OffsetVariantMap, kgl::FilteredVariantStats>
kgl::SequenceVariantFilter::getCanonicalVariants(const std::shared_ptr<const ContigDB>& contig_ptr,
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
  auto hetero_contig_ptr = modify_contig_ptr->viewFilter(HeterozygousFilter());
  size_t total_interval_variants = hetero_contig_ptr->variantCount();
  auto snp_contig_ptr = hetero_contig_ptr->viewFilter(SNPFilter());
  size_t total_snp_variants = snp_contig_ptr->variantCount();
  auto frame_shift_contig_ptr = hetero_contig_ptr->viewFilter(FrameShiftFilter());
  size_t total_frame_shift = frame_shift_contig_ptr->variantCount();

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

  FilteredVariantStats filter_stats;
  filter_stats.non_unique_count_ = non_unique_count;
  filter_stats.upstream_deleted_ = upstream_deleted;
  filter_stats.total_interval_variants_ = total_interval_variants;
  filter_stats.total_snp_variants_ = total_snp_variants;
  filter_stats.total_frame_shift_ = total_frame_shift;
  return { offset_variant_map, filter_stats };

}


void kgl::SequenceVariantFilter::canonicalVariants(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                   const OpenRightUnsigned& sequence_interval) {


  const auto [interval_map, filter_stats] = SequenceVariantFilter::getCanonicalVariants(contig_ptr, sequence_interval);
  offset_variant_map_ = std::move(interval_map);
  variant_filter_stats_ = filter_stats;

}


void kgl::SequenceVariantFilter::selectFilterType(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                  const OpenRightUnsigned& sequence_interval) {

  switch(sequenceFilterType()) {

    case SeqVariantFilterType::DEFAULT_SEQ_FILTER:
    case SeqVariantFilterType::HIGHEST_FREQ_VARIANT:
      canonicalVariants(contig_ptr, sequence_interval);
      break;

    case SeqVariantFilterType::SNP_ONLY_VARIANT:
      canonicalVariants(contig_ptr, sequence_interval);
      break;

  }

}
