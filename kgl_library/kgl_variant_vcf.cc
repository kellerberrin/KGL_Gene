//
// Created by kellerberrin on 31/08/18.
//


#include "kgl_variant_vcf.h"
#include "kgl_sequence_offset.h"
#include "kgl_variant_factory_vcf_parse_impl.h"

#include <algorithm>


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::VCFVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail and evidence()) {

    ss << evidence()->output(delimiter, output_index);

  }

  return ss.str();

}



bool kgl::VCFVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const VCFVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return genomeId() == cmp_snp->genomeId()
         and contigId() == cmp_snp->contigId()
         and phaseId() == cmp_snp->phaseId()
         and offset() == cmp_snp->offset()
         and variantType() == cmp_snp->variantType()
         and reference().getSequenceAsString() == cmp_snp->reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_snp->alternate().getSequenceAsString();

}


// Order variant types.
bool kgl::VCFVariant::lessThan(const Variant& cmp_var) const {


  if (genomeId() < cmp_var.genomeId()) {

    return true;

  } else if (contigId() < cmp_var.contigId()) {

    return true;

  } else if (contigId() > cmp_var.contigId()) {

    return false;

  } else if (phaseId() < cmp_var.phaseId()) {

    return true;

  } else if (phaseId() > cmp_var.phaseId()) {

    return false;

  } else if (offset() < cmp_var.offset()) {

    return true;

  } else if (offset() > cmp_var.offset()) {

    return false;

  } else if (variantType() < cmp_var.variantType()) {

    return true;

  } else if (variantType() > cmp_var.variantType()) {

    return false;

  }

  auto cmp_snp = dynamic_cast<const VCFVariant*>(&cmp_var);

  if (not cmp_snp) {

    // Must be a variant type == snp type.
    ExecEnv::log().error("SNPVariant::lessThan; Expected VCFVariant, got: {}", cmp_var.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  if (reference().getSequenceAsString() < cmp_snp->reference().getSequenceAsString()) {

    return true;

  } else if (reference().getSequenceAsString() > cmp_snp->reference().getSequenceAsString()) {

    return false;

  } else if (alternate().getSequenceAsString() < cmp_snp->alternate().getSequenceAsString()) {

    return true;

  }

  return false;

}


std::string kgl::VCFVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << reference().getSequenceAsString() << ">" << offsetOutput(offset(), output_index) << ">";
  ss << alternate().getSequenceAsString() << delimiter;
  ss << alternateCigar() << delimiter;

  return ss.str();

}


bool kgl::VCFVariant::mutateSequence(SignedOffset_t offset_adjust,
                                     std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                     SignedOffset_t& sequence_size_modify) const {


  sequence_size_modify = 0;

  SignedOffset_t adjusted_offset = offset() + offset_adjust;

  // Check the adjusted offset
  if (adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence_ptr->length())) {

    ExecEnv::log().error("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                         adjusted_offset, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  // Check the offset to see if variant preceded the start of the sequence.
  if (adjusted_offset < 0) {

    if (adjusted_offset + referenceSize() > 0) {

      auto ref_offset = static_cast<ContigOffset_t>(std::abs(adjusted_offset));
      auto ref_size = static_cast<ContigSize_t>(referenceSize() - ref_offset);

      if (ref_size > alternate().length()) {

        // no further action.
        return true;

      }

      std::shared_ptr<DNA5SequenceLinear> adjusted_reference = reference().unstrandedRegion(ref_offset, ref_size);

      ContigOffset_t alt_offset = size() - ref_size;

      std::shared_ptr<DNA5SequenceLinear> adjusted_alternate = alternate().unstrandedRegion(alt_offset, ref_size);

      if (not performMutation(0, dna_sequence_ptr, *adjusted_reference, *adjusted_alternate)) {

        ExecEnv::log().info("mutateSequence(), problem mutating sequence with preceeding overlapping variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true));
        ExecEnv::log().info("mutateSequence(), preceeding overlapping adj. reference: {} adj. alternate: {}",
                            adjusted_reference->getSequenceAsString(), adjusted_alternate->getSequenceAsString());

        return false;

      }

      sequence_size_modify = adjusted_alternate->length() - adjusted_reference->length();

      return true;

    } else {

      ExecEnv::log().info("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                           adjusted_offset, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

  }


  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = dna_sequence_ptr->length() - sequence_offset;

  // Check that we are not deleting beyond the end of the sequence.
  if (reference().length() < max_delete_size) {

    if (not performMutation(sequence_offset, dna_sequence_ptr, reference(), alternate())) {

      ExecEnv::log().info("mutateSequence(), problem mutating sequence; variant: {}",
                           output(' ', VariantOutputIndex::START_0_BASED, true));

      return false;

    }

    sequence_size_modify = alternate().length() - reference().length();

  } else {


    std::shared_ptr<DNA5SequenceLinear> adjusted_reference = reference().unstrandedRegion(0, max_delete_size);

    size_t alternate_length = alternateSize(max_delete_size);

    std::shared_ptr<DNA5SequenceLinear> adjusted_alternate = alternate().unstrandedRegion(0, alternate_length);

    if (not performMutation(sequence_offset, dna_sequence_ptr, *adjusted_reference, *adjusted_alternate)) {

      ExecEnv::log().info("mutateSequence(), problem mutating sequence with overlapping variant: {}",
                           output(' ', VariantOutputIndex::START_0_BASED, true));
      ExecEnv::log().info("mutateSequence(), overlapping adj. reference: {} adj. alternate: {}",
                           adjusted_reference->getSequenceAsString(), adjusted_alternate->getSequenceAsString());

      return false;

    }

    sequence_size_modify = adjusted_alternate->length() - adjusted_reference->length();

  }


  return true;

}


bool kgl::VCFVariant::performMutation(ContigOffset_t offset,
                                      std::shared_ptr<DNA5SequenceLinear> mutated_sequence_ptr,
                                      const DNA5SequenceLinear& delete_subsequence,
                                      const DNA5SequenceLinear& add_subsequence) {


  const size_t suffix_prefix = 10;
  auto reference_offset = offset;
  bool reference_check = true;

  // Check the offset
  if (offset + delete_subsequence.length() > mutated_sequence_ptr->length()) {

    ExecEnv::log().error("performMutation(), sequence offset: {} + delete sequence size: {} is out of range for mutate sequence size: {}",
                         offset, delete_subsequence.length(), mutated_sequence_ptr->length());
    return false;
  }

  for (size_t idx = 0; idx < delete_subsequence.length(); ++idx) {

    // Check the reference.
    if (delete_subsequence[idx] != mutated_sequence_ptr->at(reference_offset)) {

      ExecEnv::log().info("performMutation(), reference base: {} does not match sequence base: {} at (size: {}) offset: {}, delete (reference) offset: {}",
                          DNA5::convertToChar(delete_subsequence[idx]),
                          DNA5::convertToChar(mutated_sequence_ptr->at(reference_offset)),
                          mutated_sequence_ptr->length(),
                          offset, idx);

      SignedOffset_t prefix_offset = offset - suffix_prefix;
      ContigOffset_t p_offset = prefix_offset < 0 ? 0 : static_cast<ContigOffset_t>(prefix_offset);
      std::string prefix_string = mutated_sequence_ptr->unstrandedRegion(p_offset, ((2 * suffix_prefix)+delete_subsequence.length()))->getSequenceAsString();
      std::string sequence_string = mutated_sequence_ptr->unstrandedRegion(offset, delete_subsequence.length())->getSequenceAsString();
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
    mutated_sequence_ptr->modifyBase(offset, add_subsequence[0]);

  } else {

    // Mutate the sequence
    // Delete the reference
    if (not mutated_sequence_ptr->deleteSubSequence(offset, delete_subsequence.length())) {

      ExecEnv::log().error("performMutation(), could not delete at offset: {}, delete size: {}, sequence length: {} , variant: {}",
                           offset, delete_subsequence.length(), mutated_sequence_ptr->length());
      return false;

    }

    // Insert the alternate
    if (not mutated_sequence_ptr->insertSubSequence(offset, add_subsequence)) {

      ExecEnv::log().error("performMutation(), could not insert at offset: {}, insert size: {}, sequence length: {}",
                           offset, add_subsequence.length(), mutated_sequence_ptr->length());
      return false;
    }

  }

  return true;

}


std::string kgl::VCFVariant::alternateCigar() const {

  return ParseVCFMiscImpl::generateCigar(reference().getSequenceAsString(), alternate().getSequenceAsString());

}


size_t kgl::VCFVariant::alternateSize(size_t reference_size) const {

  CigarVector cigar_vector = ParseVCFMiscImpl::generateEditVector(reference().getSequenceAsString(), alternate().getSequenceAsString());
  return ParseVCFMiscImpl::alternateCount(reference_size, cigar_vector);

}
