//
// Created by kellerberrin on 31/08/18.
//


#include "kgl_variant_vcf.h"
#include "kgl_sequence_offset.h"


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

  ss << '\n';

  return ss.str();

}



bool kgl::VCFVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const VCFVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and phaseId() == cmp_snp->phaseId()
         and offset() == cmp_snp->offset()
         and variantType() == cmp_snp->variantType()
         and reference().getSequenceAsString() == cmp_snp->reference().getSequenceAsString()
         and alternate().getSequenceAsString() == cmp_snp->alternate().getSequenceAsString();

}


// Order variant types.
bool kgl::VCFVariant::lessThan(const Variant& cmp_var) const {


  if (contigId() < cmp_var.contigId()) {

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

  return ss.str();

}


bool kgl::VCFVariant::mutateSequence(SignedOffset_t offset_adjust,
                                     std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr,
                                     SignedOffset_t& sequence_size_modify) const {


  SignedOffset_t adjusted_offset = offset() + offset_adjust;

  // Check the offset
  if (adjusted_offset < 0 or adjusted_offset >= static_cast<SignedOffset_t>(dna_sequence_ptr->length())) {

    ExecEnv::log().error("mutateSequence(), calculated sequence offset: {} is out of range for sequence size: {}, variant: {}",
                         adjusted_offset, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
    return false;
  }

  auto sequence_offset = static_cast<ContigOffset_t>(adjusted_offset);
  ContigOffset_t max_delete_size = dna_sequence_ptr->length() - sequence_offset;
  ContigSize_t delete_size;

  // Check that we are not deleting beyond the end of the sequence.
  if (reference().length() > max_delete_size) {

    delete_size = max_delete_size;
    ExecEnv::log().info("mutateSequence(), compound deletion size: {},  offset: {}, sequence size: {}, max delete size: {}",
                         size(), sequence_offset, dna_sequence_ptr->length(), max_delete_size);

  } else {

    delete_size = reference().length();

  }

  auto reference_offset = sequence_offset;
  ContigSize_t reference_count = 0;

  for (size_t idx = 0; idx < reference().length(); ++idx) {

    if (reference_count >= delete_size) break;

    // Check the reference.
    if (reference().at(idx) != dna_sequence_ptr->at(reference_offset)) {

      ExecEnv::log().warn("mutateSequence(), reference base: {} does not match sequence base: {} at contig: {} offset: {}, reference offset: {}",
                          DNA5::convertToChar(reference().at(idx)),
                          DNA5::convertToChar(dna_sequence_ptr->at(reference_offset)),
                          contigId(), offset(), idx);

      std::string seq_reference = dna_sequence_ptr->unstrandedRegion(sequence_offset, delete_size)->getSequenceAsString();
      const ContigSize_t front_porch = 10;
      SignedOffset_t preface_offset = (sequence_offset - front_porch);
      const ContigOffset_t porch_offset = preface_offset < 0 ? 0 : (sequence_offset - front_porch);
      std::string porch_str = dna_sequence_ptr->unstrandedRegion(porch_offset, (sequence_offset - porch_offset))->getSequenceAsString();
      ExecEnv::log().warn("mutateSequence(), seq preface: {}, seq reference: {}, seq length: {}, seq index: {}, offset adjust: {}",
                          porch_str, seq_reference, dna_sequence_ptr->length(), reference_offset, offset_adjust);
      ExecEnv::log().warn("mutateSequence(), mutation variant: {}", output(' ', VariantOutputIndex::START_0_BASED, true));



    }

    ++reference_offset;
    ++reference_count;

  }

  sequence_size_modify = 0;

  // If the variant is an SNP then just modify the relevant base for performance reasons.
  if (isSNP()) {

    // Mutate the sequence
    // SNPs do not modify the sequence size.
    dna_sequence_ptr->modifyBase(sequence_offset, alternate().at(0));

  } else {

    // Mutate the sequence
    // Delete the reference
    if (not dna_sequence_ptr->deleteSubSequence(sequence_offset, delete_size)) {

      ExecEnv::log().error("mutateSequence(), could not delete at offset: {}, delete size: {}, sequence length: {} , variant: {}",
                           sequence_offset, delete_size, dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;

    }

    sequence_size_modify -= delete_size;

    // Mutate the sequence
    // Insert the alternate
    DNA5SequenceLinear insert_seq(alternate().getAlphabetString());

    if (not dna_sequence_ptr->insertSubSequence(sequence_offset, insert_seq)) {

      ExecEnv::log().error("mutateSequence(), could not insert at offset: {}, alternate size: {}, sequence length: {} , variant: {}",
                           sequence_offset, alternate().length(), dna_sequence_ptr->length(), output(' ', VariantOutputIndex::START_0_BASED, true));
      return false;
    }

    sequence_size_modify += alternate().length();

  }

  return true;

}
