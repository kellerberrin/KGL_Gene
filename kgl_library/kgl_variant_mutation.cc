//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
#include "kgl_variant_mutation.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


// Mutate coding sequences
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     const std::shared_ptr<const ContigReference>& contig_ptr,
                                     const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                                     DNA5SequenceCoding& dna_sequence) {


  // Mutate unstranded DNA and then convert to STRANDED DNA
  ContigSize_t sequence_size = coding_sequence_ptr->end() - coding_sequence_ptr->start();

  // Perform the mutation
  DNA5SequenceLinear unstranded;
  if (not mutateDNA(variant_map, contig_ptr, coding_sequence_ptr->start(), sequence_size, unstranded)) {

    ExecEnv::log().error("mutateDNA(), DNA mutation failed for contig: {}", contig_ptr->contigId());

    return false;

  }

  // Convert to stranded DNA.
  dna_sequence = kgl::SequenceOffset::mutantCodingSubSequence( coding_sequence_ptr,
                                                               unstranded,
                                                               variant_mutation_offset_,
                                                               0,
                                                               0,
                                                               coding_sequence_ptr->start());

  return true;

}

// Mutate DNA region.
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& region_variant_map,
                                     const std::shared_ptr<const ContigReference>& contig_ptr,
                                     const ContigOffset_t contig_offset,
                                     const ContigSize_t sequence_size,
                                     DNA5SequenceLinear& dna_sequence) {

  dna_sequence = contig_ptr->sequence().subSequence(contig_offset, sequence_size);

  SignedOffset_t sequence_size_modify;
  variant_mutation_offset_.clearIndelOffset();

  for (auto const& variant : region_variant_map) {

    // Adjust the mutation offset for indels.
    SignedOffset_t adjusted_offset = variant_mutation_offset_.adjustIndelOffsets(variant.second->offset());

    // Adjust the offset for the sequence offset
    adjusted_offset = adjusted_offset - contig_offset;

    // Mutate the sequence
    if (not variant.second->mutateSequence(adjusted_offset, dna_sequence, sequence_size_modify)) {

      ExecEnv::log().info("mutateDNA(), DNA mutation failed for variant: {}",
                          variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

      ExecEnv::log().info("mutateDNA(), Offset: {}, Sequence Length: {}, list of all sequence variants follows:",
                          contig_offset, dna_sequence.length());

      for(auto const& map_variant : region_variant_map) {

        ExecEnv::log().info("mutateDNA(), variant: {}",
                            map_variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

      }

    }

    // Update the mutation offset for indels.
    if (not variant_mutation_offset_.updateIndelAccounting(variant.second, sequence_size_modify)) {

      ExecEnv::log().error( "Problem updating indel mutated sequence length by variant: {}",
                            variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

    }

  }

  // Check the sequence size.
  ContigSize_t calc_sequence_size = sequence_size + variant_mutation_offset_.totalIndelOffset();

  if (calc_sequence_size != dna_sequence.length()) {

    ExecEnv::log().error("Mutated sequence length: {}, unmutated size: {}, indel adjust: {}",
                         dna_sequence.length(), sequence_size, variant_mutation_offset_.totalIndelOffset());

  }

  return true;

}


bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     const std::shared_ptr<const ContigReference>& contig_ptr,
                                     DNA5SequenceContig& contig_sequence) {

  if (not mutateDNA(variant_map, contig_ptr, 0, contig_ptr->sequence().length(), contig_sequence)) {

    ExecEnv::log().error("mutateDNA(), problem to mutating contig: {}", contig_ptr->contigId());
    return false;

  }

  return true;

}

