//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
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
    if (not variant.second->mutateSequence(adjusted_offset, dna_sequence_ptr, sequence_size_modify)) {

      ExecEnv::log().info("mutateDNA(), DNA mutation failed for variant: {}",
                          variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

      ExecEnv::log().info("mutateDNA(), Offset: {}, Sequence Length: {}, list of all sequence variants follows:",
                          sequence_offset, dna_sequence_ptr->length());

      for(auto map_variant : region_variant_map) {

        ExecEnv::log().info("mutateDNA(), sequence variant: {}",
                            map_variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

      }

    }

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


