//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
#include "kgl_variant_single.h"
#include "kgl_variant_compound.h"
#include "kgl_variant_mutation.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


// Mutate coding sequences
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     std::shared_ptr<const ContigFeatures> contig_ptr,
                                     std::shared_ptr<const CodingSequence> coding_sequence_ptr,
                                     std::shared_ptr<DNA5SequenceCoding>& dna_sequence_ptr) {

  // mutate UNSTRANDED DNA and then convert to STRANDED DNA
  // First we extract the UNSTRANDED coding region.
  std::shared_ptr<DNA5SequenceLinear>
  unstranded_ptr = contig_ptr->sequence().unstrandedRegion(coding_sequence_ptr->start(),
                                                           (coding_sequence_ptr->end() - coding_sequence_ptr->start()));

  // And mutate it.
  mutateDNA(variant_map, coding_sequence_ptr->start(), unstranded_ptr);

  // Convert to stranded DNA.
  dna_sequence_ptr = kgl::SequenceOffset::mutantCodingSubSequence(coding_sequence_ptr,
                                                                  *unstranded_ptr,
                                                                  variant_mutation_offset_,
                                                                  0,
                                                                  0,
                                                                  coding_sequence_ptr->start());

  return true;

}

// Mutate DNA region.
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& region_variant_map,
                                     ContigOffset_t sequence_offset,
                                     std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) {

  ContigSize_t sequence_size = dna_sequence_ptr->length();
  SignedOffset_t sequence_size_modify;
  variant_mutation_offset_.clearIndelOffset();

  for (auto variant : region_variant_map) {

    // Adjust the mutation offset for indels.
    SignedOffset_t adjusted_offset = variant_mutation_offset_.adjustIndelOffsets(variant.second->offset());

    // Adjust the offset for the sequence offset
    adjusted_offset = adjusted_offset - sequence_offset;

    // Mutate the sequence
    variant.second->mutateSequence(adjusted_offset, dna_sequence_ptr, sequence_size_modify);

    // Update the mutation offset for indels.
    variant_mutation_offset_.updateIndelAccounting(variant.second, sequence_size_modify);

  }

  // Check the sequence size.
  ContigSize_t calc_sequence_size = sequence_size + variant_mutation_offset_.totalIndelOffset();

  if (calc_sequence_size != dna_sequence_ptr->length()) {

    ExecEnv::log().error("Mutated sequence length: {}, unmutated size: {}, indel adjust: {}",
                         dna_sequence_ptr->length(), sequence_size, variant_mutation_offset_.totalIndelOffset());

  }

  return true;

}


// Recursively generates mutation alternatives.
// Warning - not a long function, but the coding logic is convoluted.
// Study carefully before modification.
void kgl::VariantMutation::getMutationAlternatives(std::shared_ptr<const OffsetVariantMap> variant_map_ptr,
                                                   std::vector<OffsetVariantMap>& variant_map_vector,
                                                   size_t& alternative_count,
                                                   size_t soft_limit,
                                                   size_t hard_limit) {

  alternative_count += 1;

  if (alternative_count == soft_limit) {

    ExecEnv::log().info("Soft limit of {} coding sequence mutations reached", soft_limit);

  }

  if (alternative_count > hard_limit) {

    ExecEnv::log().warn("Hard limit of {} coding sequence mutations reached, no additional mutations generated", hard_limit);
    return;

  }

  OffsetVariantMap alternative_map;
  for (auto it = variant_map_ptr->begin(); it != variant_map_ptr->end(); ++it) {

    if (not alternative_map.empty()) {

      if (alternative_map.rbegin()->second->offsetOverlap(*it->second)) {

        std::shared_ptr<OffsetVariantMap> copy_map_ptr(std::make_shared<OffsetVariantMap>(alternative_map));
        copy_map_ptr->erase(std::prev(copy_map_ptr->end()));  // Pop the last element

        for (auto copy_it = it; copy_it != variant_map_ptr->end(); ++copy_it) { // Copy the rest.

          copy_map_ptr->insert(*copy_it);

        }
        // Recursive call.
        getMutationAlternatives(copy_map_ptr, variant_map_vector, alternative_count, soft_limit, hard_limit);

      } else { // no matching offset variant.

        alternative_map.insert(*it);

      }

    } else { // empty

      alternative_map.insert(*it);

    }

  }

  variant_map_vector.push_back(alternative_map);

}

