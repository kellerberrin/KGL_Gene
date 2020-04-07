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


  // Mutate unstranded DNA and then convert to STRANDED DNA
  ContigSize_t sequence_size = coding_sequence_ptr->end() - coding_sequence_ptr->start();

  // Perform the mutation
  std::shared_ptr<DNA5SequenceLinear> unstranded_ptr;
  if (not mutateDNA(variant_map, contig_ptr, coding_sequence_ptr->start(), sequence_size, unstranded_ptr)) {

    ExecEnv::log().error("mutateDNA(), DNA mutation failed for contig: {}", contig_ptr->contigId());

    return false;

  }

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
                                     std::shared_ptr<const ContigFeatures> contig_ptr,
                                     ContigOffset_t contig_offset,
                                     ContigSize_t sequence_size,
                                     std::shared_ptr<DNA5SequenceLinear>& dna_sequence_ptr) {

  dna_sequence_ptr = contig_ptr->sequence().subSequence(contig_offset, sequence_size);

  SignedOffset_t sequence_size_modify;
  variant_mutation_offset_.clearIndelOffset();

  for (auto variant : region_variant_map) {

    // Adjust the mutation offset for indels.
    SignedOffset_t adjusted_offset = variant_mutation_offset_.adjustIndelOffsets(variant.second->offset());

    // Adjust the offset for the sequence offset
    adjusted_offset = adjusted_offset - contig_offset;

    // Mutate the sequence
    if (not variant.second->mutateSequence(adjusted_offset, dna_sequence_ptr, sequence_size_modify)) {

      ExecEnv::log().info("mutateDNA(), DNA mutation failed for variant: {}",
                          variant.second->output(' ',VariantOutputIndex::START_0_BASED, true));

      ExecEnv::log().info("mutateDNA(), Offset: {}, Sequence Length: {}, list of all sequence variants follows:",
                          contig_offset, dna_sequence_ptr->length());

      for(auto map_variant : region_variant_map) {

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

  if (calc_sequence_size != dna_sequence_ptr->length()) {

    ExecEnv::log().error("Mutated sequence length: {}, unmutated size: {}, indel adjust: {}",
                         dna_sequence_ptr->length(), sequence_size, variant_mutation_offset_.totalIndelOffset());

  }

  return true;

}


bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     std::shared_ptr<const ContigFeatures> contig_ptr,
                                     std::shared_ptr<const DNA5SequenceContig>& contig_sequence_ptr) {

  std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr;
  if (not mutateDNA(variant_map, contig_ptr, 0, contig_ptr->sequence().length(), dna_sequence_ptr)) {

    ExecEnv::log().error("mutateDNA(), problem to mutating contig: {}", contig_ptr->contigId());
    contig_sequence_ptr = contig_ptr->sequence_ptr();
    return false;

  }

  contig_sequence_ptr = std::make_shared<DNA5SequenceContig>(dna_sequence_ptr->getAlphabetString());

  return true;

}

