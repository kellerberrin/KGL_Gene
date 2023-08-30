//
// Created by kellerberrin on 2/08/23.
//

#include "kgl_mutation_region.h"


#include <memory>
#include "kgl_mutation_variant.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


// Mutate DNA region.
void kgl::ModifiedDNARegion::mutateDNA(const RegionVariantMap& region_variant_map,
                                       const std::shared_ptr<const ContigReference>& contig_ptr) {

  ContigOffset_t contig_offset = region_variant_map.variantRegion().lower();
  size_t region_size = region_variant_map.variantRegion().size();

  // Check that the contig id matches. If not, issue a warning and continue.
  if (contig_ptr->contigId() != region_variant_map.contigId()) {

    ExecEnv::log().warn("ModifiedDNARegion::mutateDNA; Supplied contig id: {} does not match variant map contig id: {}",
                        contig_ptr->contigId() , region_variant_map.contigId());

  }

  // Get the unmutated sequence.
  modified_sequence_ = contig_ptr->sequence().subSequence(contig_offset, region_size);
  // Initialize the offset structure that records indel offset changes to the original sequence.
  variant_modification_offset_.initialSequence(contig_offset, region_size);

  // For all variants modifying the sequence.
  for (auto const& [_offset, variant_ptr] : region_variant_map.variantMap()) {

    // Modifying variants must be canonical.
    if (not variant_ptr->isCanonical()) {

      ExecEnv::log().error("ModifiedDNARegion::mutateDNA; DNA mutation attempted for non-canonical variant: {}", variant_ptr->HGVS_Phase());

    }

    // Adjust the mutation offset for indels.
    SignedOffset_t adjusted_offset = variant_modification_offset_.adjustIndelOffsets(variant_ptr->offset());

    // Adjust the offset for the sequence offset
    adjusted_offset = adjusted_offset - static_cast<SignedOffset_t>(contig_offset);

    // Mutate the sequence
    auto [sequence_size_modify, result] = modifyVariant(variant_ptr, adjusted_offset);
    if (not result) {

      ExecEnv::log().warn("ModifiedDNARegion::mutateDNA; DNA mutation failed for variant: {}", variant_ptr->HGVS_Phase());
      ExecEnv::log().info("ModifiedDNARegion::mutateDNA; Offset: {}, Sequence Length: {}", contig_offset, modified_sequence_.length());

    }

    // Update the mutation offset for indels.
    if (not variant_modification_offset_.updateIndelAccounting(variant_ptr->offset(), sequence_size_modify)) {

      ExecEnv::log().error("ModifiedDNARegion::mutateDNA; Problem updating indel mutated sequence length by variant: {}", variant_ptr->HGVS_Phase());

    }

  }

  // Check the sequence size after variant mutation.
  ContigSize_t calc_sequence_size = region_size + variant_modification_offset_.totalIndelOffset();
  if (calc_sequence_size != modified_sequence_.length()) {

    ExecEnv::log().error("ModifiedDNARegion::mutateDNA; Mutated sequence length: {}, unmutated size: {}, indel adjust: {}",
                         modified_sequence_.length(), region_size, variant_modification_offset_.totalIndelOffset());

  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::pair<kgl::SignedOffset_t, bool> kgl::ModifiedDNARegion::modifyVariant(const std::shared_ptr<const Variant>& variant_ptr,
                                                                           SignedOffset_t offset_adjust) {


  SignedOffset_t adjusted_offset = static_cast<SignedOffset_t>(variant_ptr->offset()) + offset_adjust;

  // Check the adjusted offset
  if (adjusted_offset >= static_cast<SignedOffset_t>(modified_sequence_.length())) {

    ExecEnv::log().error("ModifiedDNARegion::modifyVariant; calculated sequence offset: {} is out of range for sequence size: {}",
                         adjusted_offset, modified_sequence_.length());
    return {0, false};

  }

  // Check the offset to see if variant preceded the start of the sequence.
  if (adjusted_offset < 0) {

    auto [preceeding_size_modify, result] = precedingVariant(variant_ptr, adjusted_offset);
    return { preceeding_size_modify, result};

  }

  // adjusted_offset must be positive.
  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = modified_sequence_.length() - sequence_offset;

  auto [variant_offset_modify, result] = performModification(variant_ptr, sequence_offset);
  if (not result) {

    ExecEnv::log().info("ModifiedDNARegion::modifyVariant; problem mutating sequence with variant: {}", variant_ptr->HGVS());
    return {0, false};

  }

  return {variant_offset_modify, true};

}

std::pair<kgl::SignedOffset_t, bool> kgl::ModifiedDNARegion::precedingVariant(const std::shared_ptr<const Variant>& variant_ptr,
                                                                              SignedOffset_t adjusted_offset) {

  if (adjusted_offset + variant_ptr->referenceSize() > 0) {

    auto ref_offset = static_cast<ContigOffset_t>(std::abs(adjusted_offset));
    auto ref_size = static_cast<ContigSize_t>(variant_ptr->referenceSize() - ref_offset);

    if (ref_size > variant_ptr->alternateSize()) {

      // no further action.
      return {0, true};

    }

    DNA5SequenceLinear adjusted_reference = variant_ptr->reference().subSequence(ref_offset, ref_size);

    ContigOffset_t alt_offset = variant_ptr->alternateSize() - ref_size;

    DNA5SequenceLinear adjusted_alternate = variant_ptr->alternate().subSequence(alt_offset, ref_size);

    auto [offset_adjust, result] = performModification(variant_ptr, 0);
    if (not result) {

      ExecEnv::log().warn("ModifiedDNARegion::precedingVariant; problem mutating sequence with preceeding overlapping variant");
      ExecEnv::log().info("ModifiedDNARegion::precedingVariant; preceeding overlapping adj. reference: {} adj. alternate: {}",
                          adjusted_reference.getSequenceAsString(), adjusted_alternate.getSequenceAsString());

      return {0, false};

    }

    return {offset_adjust, true};

  } else {

    ExecEnv::log().warn("ModifiedDNARegion::precedingVariant; calculated sequence offset: {} is out of range for sequence size: {}",
                        adjusted_offset, modified_sequence_.length());
    return {0, false};

  }

}


std::pair<kgl::SignedOffset_t, bool> kgl::ModifiedDNARegion::performModification(const std::shared_ptr<const Variant>& variant_ptr,
                                                                                 ContigOffset_t variant_offset) {


  SignedOffset_t sequence_size_modify{0};
  const size_t suffix_prefix{10};
  auto reference_offset{variant_offset};
  bool reference_check{true};

  // Check the canonical_offset
  if (variant_offset + variant_ptr->referenceSize() > modified_sequence_.length()) {

    ExecEnv::log().error("ModifiedDNARegion::performModification; sequence canonical_offset: {} + delete sequence size: {} is out of range for mutate sequence size: {}",
                         variant_offset, variant_ptr->referenceSize(), modified_sequence_.length());
    return {0, false};
  }

  for (size_t idx = 0; idx < variant_ptr->referenceSize(); ++idx) {

    // Check the reference.
    if (variant_ptr->reference().at(idx) != modified_sequence_.at(reference_offset)) {

      ExecEnv::log().info("ModifiedDNARegion::performModification; reference base: {} does not match sequence base: {} at (size: {}) canonical_offset: {}, delete (reference) canonical_offset: {}",
                          DNA5::convertToChar(variant_ptr->reference().at(idx)),
                          DNA5::convertToChar(modified_sequence_.at(reference_offset)),
                          modified_sequence_.length(),
                          variant_offset, idx);

      SignedOffset_t prefix_offset = variant_offset - suffix_prefix;
      ContigOffset_t p_offset = prefix_offset < 0 ? 0 : static_cast<ContigOffset_t>(prefix_offset);
      std::string prefix_string = modified_sequence_.subSequence(p_offset,
                                                                 ((2 * suffix_prefix) + variant_ptr->referenceSize())).getSequenceAsString();
      std::string sequence_string = modified_sequence_.subSequence(variant_offset, variant_ptr->referenceSize()).getSequenceAsString();
      ExecEnv::log().info("ModifiedDNARegion::performModification; seq (-10/ref/+10) section: {}, seq: {}, reference: {}, alternate: {}",
                          prefix_string, sequence_string, variant_ptr->reference().getSequenceAsString(), variant_ptr->alternate().getSequenceAsString());

      reference_check =false;

    }

    ++reference_offset;

  }

  if (not reference_check) return {0, false};

  if (variant_ptr->isSNP()) {

    // Mutate the sequence
    // SNPs do not modify the sequence size.
    sequence_size_modify = 0;
    if (not modified_sequence_.modifyBase(variant_offset, variant_ptr->alternate().at(0))) {

      ExecEnv::log().error("ModifiedDNARegion::performModification; could not modify base (SNP) at canonical_offset: {}", variant_offset);
      return {0, false};

    }

  } else {

    // Mutate the sequence
    // Delete the reference
    sequence_size_modify = sequence_size_modify - variant_ptr->referenceSize();
    if (not modified_sequence_.deleteSubSequence(variant_offset, variant_ptr->referenceSize())) {

      ExecEnv::log().error("ModifiedDNARegion::performModification; could not delete at canonical_offset: {}, delete size: {}, sequence length: {} , variant: {}",
                           variant_offset, variant_ptr->referenceSize(), modified_sequence_.length());
      return {0, false};

    }

    // Insert the alternate
    sequence_size_modify = sequence_size_modify + variant_ptr->alternateSize();
    if (not modified_sequence_.insertSubSequence(variant_offset, variant_ptr->alternate())) {

      ExecEnv::log().error("ModifiedDNARegion::performModification; could not insert at canonical_offset: {}, insert size: {}, sequence length: {}",
                           variant_offset, variant_ptr->alternateSize(), modified_sequence_.length());
      return {0, false};
    }

  }

  return {sequence_size_modify, true};

}
