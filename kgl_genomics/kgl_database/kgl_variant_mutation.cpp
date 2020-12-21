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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::Variant::mutateSequence( SignedOffset_t offset_adjust,
                                   DNA5SequenceLinear& dna_sequence,
                                   SignedOffset_t& sequence_size_modify) const {


  sequence_size_modify = 0;

  SignedOffset_t adjusted_offset = offset() + offset_adjust;

  // Check the adjusted offset
  if (adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence.length())) {

    ExecEnv::log().error("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                         adjusted_offset, dna_sequence.length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  // Check the offset to see if variant preceded the start of the sequence.
  if (adjusted_offset < 0) {

    return preceedingMutation(adjusted_offset, dna_sequence, sequence_size_modify);

  }

  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = dna_sequence.length() - sequence_offset;

  // Check that we are not deleting beyond the end of the sequence.
  if (reference().length() < max_delete_size) {

    if (not performMutation(sequence_offset, dna_sequence, reference(), alternate())) {

      ExecEnv::log().info("mutateSequence(), problem mutating sequence; variant: {}",
                          output(' ', VariantOutputIndex::START_0_BASED, true));

      return false;

    }

    sequence_size_modify = alternate().length() - reference().length();

  } else {


    DNA5SequenceLinear adjusted_reference = reference().subSequence(0, max_delete_size);

    size_t alternate_length = alternateSize(max_delete_size);

    DNA5SequenceLinear adjusted_alternate = alternate().subSequence(0, alternate_length);

    if (not performMutation(sequence_offset, dna_sequence, adjusted_reference, adjusted_alternate)) {

      ExecEnv::log().info("mutateSequence(), problem mutating sequence with overlapping variant: {}",
                          output(' ', VariantOutputIndex::START_0_BASED, true));
      ExecEnv::log().info("mutateSequence(), overlapping adj. reference: {} adj. alternate: {}",
                          adjusted_reference.getSequenceAsString(), adjusted_alternate.getSequenceAsString());

      return false;

    }

    sequence_size_modify = adjusted_alternate.length() - adjusted_reference.length();

  }


  return true;

}


bool kgl::Variant::preceedingMutation( SignedOffset_t adjusted_offset,
                                       DNA5SequenceLinear& dna_sequence,
                                       SignedOffset_t& sequence_size_modify) const {

  if (adjusted_offset + referenceSize() > 0) {

    auto ref_offset = static_cast<ContigOffset_t>(std::abs(adjusted_offset));
    auto ref_size = static_cast<ContigSize_t>(referenceSize() - ref_offset);

    if (ref_size > alternate().length()) {

      // no further action.
      return true;

    }

    DNA5SequenceLinear adjusted_reference = reference().subSequence(ref_offset, ref_size);

    ContigOffset_t alt_offset = alternateSize() - ref_size;

    DNA5SequenceLinear adjusted_alternate = alternate().subSequence(alt_offset, ref_size);

    if (not performMutation(0, dna_sequence, adjusted_reference, adjusted_alternate)) {

      ExecEnv::log().info("mutateSequence(), problem mutating sequence with preceeding overlapping variant: {}",
                          output(' ', VariantOutputIndex::START_0_BASED, true));
      ExecEnv::log().info("mutateSequence(), preceeding overlapping adj. reference: {} adj. alternate: {}",
                          adjusted_reference.getSequenceAsString(), adjusted_alternate.getSequenceAsString());

      return false;

    }

    sequence_size_modify = adjusted_alternate.length() - adjusted_reference.length();

    return true;

  } else {

    ExecEnv::log().info("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                        adjusted_offset, dna_sequence.length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;

  }

}



bool kgl::Variant::performMutation( ContigOffset_t offset,
                                    DNA5SequenceLinear& mutated_sequence,
                                    const DNA5SequenceLinear& delete_subsequence,
                                    const DNA5SequenceLinear& add_subsequence) {


  const size_t suffix_prefix = 10;
  auto reference_offset = offset;
  bool reference_check = true;

  // Check the offset
  if (offset + delete_subsequence.length() > mutated_sequence.length()) {

    ExecEnv::log().error("performMutation(), sequence offset: {} + delete sequence size: {} is out of range for mutate sequence size: {}",
                         offset, delete_subsequence.length(), mutated_sequence.length());
    return false;
  }

  for (size_t idx = 0; idx < delete_subsequence.length(); ++idx) {

    // Check the reference.
    if (delete_subsequence[idx] != mutated_sequence.at(reference_offset)) {

      ExecEnv::log().info("performMutation(), reference base: {} does not match sequence base: {} at (size: {}) offset: {}, delete (reference) offset: {}",
                          DNA5::convertToChar(delete_subsequence[idx]),
                          DNA5::convertToChar(mutated_sequence.at(reference_offset)),
                          mutated_sequence.length(),
                          offset, idx);

      SignedOffset_t prefix_offset = offset - suffix_prefix;
      ContigOffset_t p_offset = prefix_offset < 0 ? 0 : static_cast<ContigOffset_t>(prefix_offset);
      std::string prefix_string = mutated_sequence.subSequence(p_offset,
                                                               ((2 * suffix_prefix) + delete_subsequence.length())).getSequenceAsString();
      std::string sequence_string = mutated_sequence.subSequence(offset, delete_subsequence.length()).getSequenceAsString();
      ExecEnv::log().info("performMutation(), seq (-10/ref/+10) section: {}, seq: {}, reference: {}, alternate: {}",
                          prefix_string, sequence_string, delete_subsequence.getSequenceAsString(), add_subsequence.getSequenceAsString());

      reference_check =false;

    }

    ++reference_offset;

  }

  if (not reference_check) return false;

  if (delete_subsequence.length() == 1 and add_subsequence.length() == 1) {

    // Mutate the sequence
    // SNPs do not modify the sequence size.
    if (not mutated_sequence.modifyBase(offset, add_subsequence[0])) {

      ExecEnv::log().error("performMutation(), could not modify base (SNP) at offset: {}", offset);
      return false;

    }

  } else {

    // Mutate the sequence
    // Delete the reference
    if (not mutated_sequence.deleteSubSequence(offset, delete_subsequence.length())) {

      ExecEnv::log().error("performMutation(), could not delete at offset: {}, delete size: {}, sequence length: {} , variant: {}",
                           offset, delete_subsequence.length(), mutated_sequence.length());
      return false;

    }

    // Insert the alternate
    if (not mutated_sequence.insertSubSequence(offset, add_subsequence)) {

      ExecEnv::log().error("performMutation(), could not insert at offset: {}, insert size: {}, sequence length: {}",
                           offset, add_subsequence.length(), mutated_sequence.length());
      return false;
    }

  }

  return true;

}
