//
// Created by kellerberrin on 26/07/23.
//

#include "kgl_mutation_offset.h"


namespace kgl = kellerberrin::genome;



void kgl::MutationOffset::checkUpstreamDeletion(OffsetVariantMap& variant_map) {

  for (auto iter = variant_map.begin(); iter != variant_map.end(); ++iter) {

    // Ignore the first entry.
    if (iter == variant_map.begin()) {

      continue;

    }

    auto const& [previous_offset, previous_variant_ptr] = *std::ranges::prev(iter, 1, variant_map.begin());
    auto const& [offset, variant_ptr] = *iter;

    int64_t delete_size = previous_variant_ptr->referenceSize() - previous_variant_ptr->alternateSize();
    int64_t offset_gap = variant_ptr->offset() - previous_variant_ptr->offset();

    if (delete_size >= offset_gap) {

      //      ExecEnv::log().info("ContigDB::checkUpstreamDeletion(), Upstream deletion detected: {}",
      //                           std::prev(iter)->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      //      ExecEnv::log().info("ContigDB::checkUpstreamDeletion(), Downstream variant removed from mutation: {}",
      //                           iter->second->output(' ', VariantOutputIndex::START_0_BASED, true));

      iter = variant_map.erase(iter);

      iter = std::ranges::prev(iter, 1, variant_map.begin()); // reset back to the previous iterator.

    }

  }

}



bool kgl::MutationOffset::getSortedVariants( const std::shared_ptr<const ContigDB>& contig_ptr,
                                             VariantPhase phase,
                                             ContigOffset_t start,
                                             ContigOffset_t end,
                                             OffsetVariantMap& variant_map) {


  auto lower_bound = contig_ptr->getMap().lower_bound(start);
  auto upper_bound = contig_ptr->getMap().upper_bound(end-1); //  [start, end)

  // If there is a prior variant that overlaps the start address, then push this onto the variant map.
  if (lower_bound != contig_ptr->getMap().end() and lower_bound != contig_ptr->getMap().begin()) {

    auto previous_offset_ptr = std::ranges::prev(lower_bound, 1, contig_ptr->getMap().begin());

    auto const& [offset, offset_db_ptr] = *previous_offset_ptr;

    const OffsetDBArray& previous_offset_variants = offset_db_ptr->getVariantArray();

    for (const auto& variant_ptr : previous_offset_variants) {

      if (variant_ptr->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        if (variant_ptr->offset() + variant_ptr->referenceSize() > start) {

          variant_map.emplace(variant_ptr->offset(), variant_ptr);

        }

      }

    }

  }

  // Get all variants between the lower and upper bounds.
  for (auto it = lower_bound; it != upper_bound; ++it) {

    auto const& [offset, offset_db_ptr] = *it;

    const OffsetDBArray& previous_offset_variants = offset_db_ptr->getVariantArray();

    for (const auto& variant_ptr : previous_offset_variants) {

      if (variant_ptr->phaseId() == phase or phase == VariantPhase::UNPHASED) {

        variant_map.emplace(offset, variant_ptr);

      }

    }

  }

  checkUpstreamDeletion(variant_map);

  return true;

}




bool kgl::MutationOffset::getSortedVariants( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                             ContigId_t contig_id,
                                             VariantPhase phase,
                                             ContigOffset_t start,
                                             ContigOffset_t end,
                                             OffsetVariantMap &variant_map) {


  auto find_iter = genome_ptr->getMap().find(contig_id);

  if (find_iter == genome_ptr->getMap().end()) {

    ExecEnv::log().error("MutationOffset::getSortedVariants; Contig Id: {} not found in Genome Variant: {}", contig_id, genome_ptr->genomeId());
    return false;

  }

  std::shared_ptr<ContigDB> contig_ptr = find_iter->second;

  return MutationOffset::getSortedVariants(contig_ptr, phase, start, end, variant_map);

}


