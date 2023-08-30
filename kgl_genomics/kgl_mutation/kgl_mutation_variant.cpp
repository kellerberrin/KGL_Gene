//
// Created by kellerberrin on 22/12/17.
//

#include <memory>
#include "kgl_mutation_variant.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


// Mutate coding sequences
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& variant_map,
                                     const std::shared_ptr<const ContigReference>& contig_ptr,
                                     const std::shared_ptr<const TranscriptionSequence>& coding_sequence_ptr,
                                     DNA5SequenceCoding& dna_sequence) {


  // Mutate mutated_unstranded_sequence DNA and then convert to STRANDED DNA
  ContigSize_t sequence_size = coding_sequence_ptr->end() - coding_sequence_ptr->start();

  // Perform the mutation
  DNA5SequenceLinear mutated_unstranded_sequence;
  if (not mutateDNA(variant_map, contig_ptr, coding_sequence_ptr->start(), sequence_size, mutated_unstranded_sequence)) {

    ExecEnv::log().error("VariantMutation::mutateDNA; DNA mutation failed for contig: {}", contig_ptr->contigId());

    return false;

  }

  // Convert to stranded DNA.
  dna_sequence = kgl::SequenceOffset::mutantCodingSubSequence(coding_sequence_ptr,
                                                              mutated_unstranded_sequence,
                                                              variant_mutation_offset_,
                                                              0,
                                                              0,
                                                              coding_sequence_ptr->start());

  return true;

}

// Mutate DNA region.
bool kgl::VariantMutation::mutateDNA(const OffsetVariantMap& region_variant_map,
                                     const std::shared_ptr<const ContigReference>& contig_ptr,
                                     ContigOffset_t contig_offset,
                                     ContigSize_t sequence_size,
                                     DNA5SequenceLinear& dna_sequence) {

  // Get the unmutated sequence.
  dna_sequence = contig_ptr->sequence().subSequence(contig_offset, sequence_size);

  // Clear the indel variant offset map and set the sequence offset and size.
  variant_mutation_offset_.initialSequence(contig_offset, sequence_size);

  // For all variants modifying the sequence.
  for (auto const& [map_offset, variant_ptr] : region_variant_map) {

    // Must be canonical; variants in the OffsetVariantMap should be pre-converted to canonical format.
    if (not variant_ptr->isCanonical()) {

      ExecEnv::log().error("VariantMutation::mutateDNA; DNA mutation attempted for non-canonical variant: {}", variant_ptr->HGVS_Phase());
      return false;

    }

    // Adjust the mutation offset for indels.
    SignedOffset_t adjusted_offset = variant_mutation_offset_.adjustIndelOffsets(variant_ptr->offset());

    // Adjust the offset for the sequence offset. This gives the adjusted from the beginning of the dna_sequence.
    adjusted_offset = adjusted_offset - static_cast<SignedOffset_t>(contig_offset);

    // Mutate the sequence
    SignedOffset_t sequence_size_modify{0};
    if (not variantMutate(variant_ptr,
                          adjusted_offset,
                          dna_sequence,
                          sequence_size_modify)) {

      ExecEnv::log().warn("VariantMutation::mutateDNA; DNA mutation failed for variant: {}", variant_ptr->HGVS_Phase());
      ExecEnv::log().info("VariantMutation::mutateDNA; Offset: {}, Sequence Length: {}", contig_offset, dna_sequence.length());

    }

    // Update the mutation offset for indels.
    if (not variant_mutation_offset_.updateIndelAccounting(variant_ptr->offset(), sequence_size_modify)) {

      ExecEnv::log().error("VariantMutation::mutateDNA; Problem updating indel mutated sequence length by variant: {}", variant_ptr->HGVS_Phase());

    }

  }

  // Check the sequence size.
  ContigSize_t calc_sequence_size = sequence_size + variant_mutation_offset_.totalIndelOffset();

  if (calc_sequence_size != dna_sequence.length()) {

    ExecEnv::log().error("VariantMutation::mutateDNA; Mutated sequence length: {}, unmutated size: {}, indel adjust: {}",
                         dna_sequence.length(), sequence_size, variant_mutation_offset_.totalIndelOffset());

  }

  return true;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VariantMutation::variantMutate( const std::shared_ptr<const Variant>& variant_ptr,
                                         SignedOffset_t offset_adjust,
                                         DNA5SequenceLinear& dna_sequence,
                                         SignedOffset_t& sequence_size_modify) {

  sequence_size_modify = 0;

  // Adjusted offset is relative to the dna_sequence.
  SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(variant_ptr->offset()) + offset_adjust;

  // Check the adjusted offset. An adjusted sequence beyond the end of the dna_sequence is an error.
  if (adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence.length())) {

    ExecEnv::log().error("VariantMutation::variantMutate; calculated sequence offset: {} is out of range for sequence size: {}",adjusted_offset, dna_sequence.length());
    return false;

  }

  // Check the offset to see if variant preceded the start of the sequence.
  if (adjusted_offset < 0) {

    return preceedingMutation(variant_ptr, adjusted_offset, dna_sequence, sequence_size_modify);

  }

  // adjusted_offset must be positive.
  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = dna_sequence.length() - sequence_offset;

  // Check that we are not deleting beyond the end of the sequence.
  if (variant_ptr->referenceSize() < max_delete_size) {

    if (not performMutation(sequence_offset, dna_sequence, variant_ptr->reference(), variant_ptr->alternate(), sequence_size_modify)) {

      ExecEnv::log().info("VariantMutation::variantMutate; problem mutating sequence");
      return false;

    }

  } else {

    if (max_delete_size > 1 and variant_ptr->alternateSize() == 1) {

      DNA5SequenceLinear adjusted_reference = variant_ptr->reference().subSequence(0, max_delete_size);

      if (not performMutation(sequence_offset, dna_sequence, adjusted_reference, variant_ptr->alternate(), sequence_size_modify)) {

        ExecEnv::log().warn("VariantMutation::variantMutate; problem mutating sequence with overlapping variant");
        ExecEnv::log().info("VariantMutation::variantMutate; overlapping adj. reference: {},  alternate: {}",
                            adjusted_reference.getSequenceAsString(),
                            variant_ptr->alternate().getSequenceAsString());
        return false;

      }

    } else {

      ExecEnv::log().warn("VariantMutation::variantMutate; max delete size: {},  overlapping reference: {} alternate: {}",
                          max_delete_size,
                          variant_ptr->reference().getSequenceAsString(),
                          variant_ptr->alternate().getSequenceAsString());
      return false;

    }

  }

  return true;

}

// The adjusted offset precedes the dna_sequence.
bool kgl::VariantMutation::preceedingMutation( const std::shared_ptr<const Variant>& variant_ptr,
                                               SignedOffset_t adjusted_offset,
                                               DNA5SequenceLinear& dna_sequence,
                                               SignedOffset_t& sequence_size_modify) {

  if (adjusted_offset + variant_ptr->referenceSize() > 0) {

    auto ref_offset = static_cast<ContigOffset_t>(std::abs(adjusted_offset));
    auto ref_size = static_cast<ContigSize_t>(variant_ptr->referenceSize() - ref_offset);

    if (ref_size > variant_ptr->alternateSize()) {

      // no further action.
      return true;

    }

    DNA5SequenceLinear adjusted_reference = variant_ptr->reference().subSequence(ref_offset, ref_size);

    ContigOffset_t alt_offset = variant_ptr->alternateSize() - ref_size;

    DNA5SequenceLinear adjusted_alternate = variant_ptr->alternate().subSequence(alt_offset, ref_size);

    if (not performMutation(0, dna_sequence, adjusted_reference, adjusted_alternate, sequence_size_modify)) {

      ExecEnv::log().warn("VariantMutation::variantMutate; problem mutating sequence with preceeding overlapping variant");
      ExecEnv::log().info("VariantMutation::variantMutate; preceeding overlapping adj. reference: {} adj. alternate: {}",
                          adjusted_reference.getSequenceAsString(), adjusted_alternate.getSequenceAsString());

      return false;

    }

    return true;

  } else {

    ExecEnv::log().warn("VariantMutation::variantMutate; calculated sequence offset: {} is out of range for sequence size: {}", adjusted_offset, dna_sequence.length());
    return false;

  }

}


bool kgl::VariantMutation::performMutation( ContigOffset_t variant_offset,
                                            DNA5SequenceLinear& mutated_sequence,
                                            const DNA5SequenceLinear& delete_subsequence,
                                            const DNA5SequenceLinear& add_subsequence,
                                            SignedOffset_t& sequence_size_modify) {


  const size_t suffix_prefix{10};
  auto reference_offset{variant_offset};
  bool reference_check{true};

  // Check the canonical_offset
  if (variant_offset + delete_subsequence.length() > mutated_sequence.length()) {

    ExecEnv::log().error("VariantMutation::performModification; sequence canonical_offset: {} + delete sequence size: {} is out of range for mutate sequence size: {}",
                         variant_offset, delete_subsequence.length(), mutated_sequence.length());
    return false;
  }

  for (size_t idx = 0; idx < delete_subsequence.length(); ++idx) {

    // Check the reference.
    if (delete_subsequence[idx] != mutated_sequence.at(reference_offset)) {

      ExecEnv::log().info("VariantMutation::performModification; reference base: {} does not match sequence base: {} at (size: {}) canonical_offset: {}, delete (reference) canonical_offset: {}",
                          DNA5::convertToChar(delete_subsequence[idx]),
                          DNA5::convertToChar(mutated_sequence.at(reference_offset)),
                          mutated_sequence.length(),
                          variant_offset, idx);

      SignedOffset_t prefix_offset = static_cast<SignedOffset_t>(variant_offset) - suffix_prefix;
      ContigOffset_t p_offset = prefix_offset < 0 ? 0 : static_cast<ContigOffset_t>(prefix_offset);
      std::string prefix_string = mutated_sequence.subSequence(p_offset,
                                                               ((2 * suffix_prefix) + delete_subsequence.length())).getSequenceAsString();
      std::string sequence_string = mutated_sequence.subSequence(variant_offset, delete_subsequence.length()).getSequenceAsString();
      ExecEnv::log().info("VariantMutation::performModification; seq (-10/ref/+10) section: {}, seq: {}, reference: {}, alternate: {}",
                          prefix_string, sequence_string, delete_subsequence.getSequenceAsString(), add_subsequence.getSequenceAsString());

      reference_check =false;

    }

    ++reference_offset;

  }

  if (not reference_check) return false;

  if (delete_subsequence.length() == 1 and add_subsequence.length() == 1) {

    // Mutate the sequence
    // SNPs do not modify the sequence size.
    sequence_size_modify = 0;
    if (not mutated_sequence.modifyBase(variant_offset, add_subsequence[0])) {

      ExecEnv::log().error("VariantMutation::performModification; could not modify base (SNP) at canonical_offset: {}", variant_offset);
      return false;

    }

  } else {

    // Mutate the sequence
    // Delete the reference
    sequence_size_modify = sequence_size_modify - delete_subsequence.length();
    if (not mutated_sequence.deleteSubSequence(variant_offset, delete_subsequence.length())) {

      ExecEnv::log().error("VariantMutation::performModification; could not delete at canonical_offset: {}, delete size: {}, sequence length: {} , variant: {}",
                           variant_offset, delete_subsequence.length(), mutated_sequence.length());
      return false;

    }

    // Insert the alternate
    sequence_size_modify = sequence_size_modify + add_subsequence.length();
    if (not mutated_sequence.insertSubSequence(variant_offset, add_subsequence)) {

      ExecEnv::log().error("VariantMutation::performModification; could not insert at canonical_offset: {}, insert size: {}, sequence length: {}",
                           variant_offset, add_subsequence.length(), mutated_sequence.length());
      return false;
    }

  }

  return true;

}
