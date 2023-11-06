//
// Created by kellerberrin on 14/08/18.
//

#include "kgl_distance_genome.h"


namespace kgl = kellerberrin::genome;

// This routine counts the the number of variants that are in one genome but not the other - the variant difference.
// Unfortunately the logic here is somewhat complex, be careful when modifying.
kgl::DistanceType_t kgl::GenomeDistance::distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const {

  std::shared_ptr<const GenomeDistance> node_ptr = std::dynamic_pointer_cast<const GenomeDistance>(distance_node);
  if (not node_ptr) {

    ExecEnv::log().error("Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  ExecEnv::log().info("Counting different variants between genome: {} and genome: {}",
                      genome_db_ptr_->genomeId(), node_ptr->genome_db_ptr_->genomeId());

  size_t distance = 0;

  // iterate through all the contigs and count mismatching variants.
  for (auto contig : genome_db_ptr_->getMap()) {

    auto cmp_contig_iter = node_ptr->genome_db_ptr_->getMap().find(contig.first);

    // Corresponding contig_ref_ptr not found.
    if (cmp_contig_iter == node_ptr->genome_db_ptr_->getMap().end()) {

      distance += contig.second->variantCount();
      continue;

    }

    // Set to first offset
    auto cmp_offset_iter = cmp_contig_iter->second->getMap().begin();

    // Scroll through offsets
    for(const auto& offset : contig.second->getMap()) {

      // If at end of cmp offsets then add the offset size.
      if (cmp_offset_iter == cmp_contig_iter->second->getMap().end()) {

        distance += offset.second->getVariantArray().size();

      } else if (offset.first == cmp_offset_iter->first) {
      // If the offsets are equal then calculate how many variants match

        size_t variants_found = 0;
        for (auto variant : offset.second->getVariantArray()) {

          for (auto cmp_variant : cmp_offset_iter->second->getVariantArray()) {

            if (variant->equivalent(*(cmp_variant))) {

              ++variants_found;
              break;

            }

          }

        }

        size_t max_size = offset.second->getVariantArray().size() >= cmp_offset_iter->second->getVariantArray().size()
           ? offset.second->getVariantArray().size() : cmp_offset_iter->second->getVariantArray().size();
        distance += (max_size - variants_found);
        // increment to next cmp offset.
        ++cmp_offset_iter;

      } else if (offset.first > cmp_offset_iter->first) {
      // If offset is greater than cmp_offset, then add all variants and increment the cmp offset.
      // Until the cmp_offset is equal to or greater than the offset.
        while (offset.first > cmp_offset_iter->first) {

          distance += cmp_offset_iter->second->getVariantArray().size();

          // increment to next cmp offset.
          ++cmp_offset_iter;

          // Check for the end condition.
          if (cmp_offset_iter == cmp_contig_iter->second->getMap().end()) {

            break;

          }

          if (offset.first == cmp_offset_iter->first) {

            // If the offsets are equal then calculate how many variants match
            size_t variants_found = 0;
            for (auto variant : offset.second->getVariantArray()) {

              for (auto cmp_variant : cmp_offset_iter->second->getVariantArray()) {

                if (variant->equivalent(*(cmp_variant))) {

                  ++variants_found;
                  break;

                }

              }

            }

            size_t max_size = offset.second->getVariantArray().size() >= cmp_offset_iter->second->getVariantArray().size()
              ? offset.second->getVariantArray().size() : cmp_offset_iter->second->getVariantArray().size();
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
        distance += offset.second->getVariantArray().size();

      }

    } // for all offsets.

    // If any cmp offsets remain then add these to the distance
    while (cmp_offset_iter != cmp_contig_iter->second->getMap().end()) {

      distance += cmp_offset_iter->second->getVariantArray().size();
      ++cmp_offset_iter;

    }

  }

  ExecEnv::log().info("GenomeDistance; counted {} different variants between genome: {} and genome: {}",
                      distance, genome_db_ptr_->genomeId(), node_ptr->genome_db_ptr_->genomeId());

  return static_cast<DistanceType_t >(distance);

}
