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

// Returns a map of unique canonical variants
// Also returns the number of multiple variants found at each offset which are filtered to a single variant.
std::pair<kgl::OffsetVariantMap, size_t> kgl::MutationOffset::getCanonicalVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                                  ContigOffset_t start,
                                                                                  ContigOffset_t end) {

  OffsetVariantMap offset_variant_map;

  // Filter the variants plus a margin of 200 nucleotides at the beginning of the specified region.
  // To allow for any offset() change when the variants are converted to canonical.
  ContigOffset_t lower = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(start)-NUCLEOTIDE_CANONICAL_MARGIN);
  // Filter to just variants in the region plus offset margins.
  auto region_contig_ptr = contig_ptr->viewFilter(ContigRegionFilter(lower, end));
  // Remove any duplicate (homozygous) identical variants.
  auto no_homo_contig_ptr = region_contig_ptr->viewFilter(UniqueUnphasedFilter());
  // Convert to canonical variants.
  auto canonical_contig_ptr = no_homo_contig_ptr->canonicalContig();
  // Filter to just the variants that will modify the specified region [start, end).
  auto modify_contig_ptr = canonical_contig_ptr->viewFilter(ContigModifyFilter(start, end));
  // Get the count of variants modifying the region [start, end).
  size_t modify_count = modify_contig_ptr->variantCount();
  // Remove multiple variants, that is minor alleles (SNP) that are not homozygous that share the same offset.
  // Note that if an indel occurs at the same offset as an SNP, this is not a problem, since in canonical form ('1MnI' or '1MnD'),
  // the indel actually modifies the next (+1) offset.
  auto unique_contig_ptr = modify_contig_ptr->viewFilter(FrequencyUniqueFilter());

  // Finally move the unique modifying variants to the offset map.
  for (auto const& [offset, offset_ptr] : unique_contig_ptr->getMap()) {

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

        ExecEnv::log().error("MutationOffset::getCanonicalVariants; duplicate variant found: {}", variant_ptr->HGVS());

      }

    } // For all variants at each offset.

  } // for all offsets.

  size_t non_unique_count = modify_count - offset_variant_map.size();

  return { offset_variant_map, non_unique_count };

}

