//
// Created by kellerberrin on 14/08/18.
//

#include "kgl_upgma_unphased.h"


namespace kgl = kellerberrin::genome;

// This routine counts the the number of variants that are in one genome but not the other - the variant difference.
// Unfortunately the logic here is somewhat complex, be careful when modifying.
kgl::DistanceType_t kgl::UPGMAUnphasedDistance::distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const {

  std::shared_ptr<const UPGMAUnphasedDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAUnphasedDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("UPGMAUnphasedDistance::distance(), Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  ExecEnv::log().info("UPGMAUnphasedDistance; counting different variants between genome: {} and genome: {}",
                      genome_variant_ptr_->genomeId(), node_ptr->genome_variant_ptr_->genomeId());

  size_t distance = 0;

  // iterate through all the contigs and count mismatching variants.
  for (auto contig : genome_variant_ptr_->getMap()) {

    auto cmp_contig_iter = node_ptr->genome_variant_ptr_->getMap().find(contig.first);

    // Corresponding contig not found.
    if (cmp_contig_iter == node_ptr->genome_variant_ptr_->getMap().end()) {

      distance += contig.second->variantCount();
      continue;

    }

    // Set to first offset
    auto cmp_offset_iter = cmp_contig_iter->second->getMap().begin();

    // Scroll through offsets
    for(auto offset : contig.second->getMap()) {

      // If at end of cmp offsets then add the offset size.
      if (cmp_offset_iter == cmp_contig_iter->second->getMap().end()) {

        distance += offset.second.size();

      } else if (offset.first == cmp_offset_iter->first) {
      // If the offsets are equal then calculate how many variants match

        size_t variants_found = 0;
        for (auto variant : offset.second) {

          for (auto cmp_variant : cmp_offset_iter->second) {

            if (variant.first->equivalent(*(cmp_variant.first))) {

              ++variants_found;
              break;

            }

          }

        }

        size_t max_size = offset.second.size() >= cmp_offset_iter->second.size() ? offset.second.size() : cmp_offset_iter->second.size();
        distance += (max_size - variants_found);
        // increment to next cmp offset.
        ++cmp_offset_iter;

      } else if (offset.first > cmp_offset_iter->first) {
      // If offset is greater than cmp_offset, then add all variants and increment the cmp offset.
      // Until the cmp_offset is equal to or greater than the offset.
        while (offset.first > cmp_offset_iter->first) {

          distance += cmp_offset_iter->second.size();

          // increment to next cmp offset.
          ++cmp_offset_iter;

          // Check for the end condition.
          if (cmp_offset_iter == cmp_contig_iter->second->getMap().end()) {

            break;

          }

          if (offset.first == cmp_offset_iter->first) {

            // If the offsets are equal then calculate how many variants match
            size_t variants_found = 0;
            for (auto variant : offset.second) {

              for (auto cmp_variant : cmp_offset_iter->second) {

                if (variant.first->equivalent(*(cmp_variant.first))) {

                  ++variants_found;
                  break;

                }

              }

            }

            size_t max_size = offset.second.size() >= cmp_offset_iter->second.size() ? offset.second.size() : cmp_offset_iter->second.size();
            distance += (max_size - variants_found);
            // increment to next cmp offset.
            ++cmp_offset_iter;

          }
          // Check for the end condition again.
          if (cmp_offset_iter == cmp_contig_iter->second->getMap().end()) {

            break;

          }

        }

      } else if (offset.first < cmp_offset_iter->first) {
      // If the offset is less than the cmp offset, then add all variants.
        distance += offset.second.size();

      }

    } // for all offsets.

    // If any cmp offsets remain then add these to the distance
    while (cmp_offset_iter != cmp_contig_iter->second->getMap().end()) {

      distance += cmp_offset_iter->second.size();
      ++cmp_offset_iter;

    }

  }

  ExecEnv::log().info("UPGMAUnphasedDistance; counted {} different variants between genome: {} and genome: {}",
                      distance, genome_variant_ptr_->genomeId(), node_ptr->genome_variant_ptr_->genomeId());

  return static_cast<DistanceType_t >(distance);

}
