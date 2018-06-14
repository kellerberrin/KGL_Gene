//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_single.h"
#include "kgl_sequence_offset.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::string kgl::SNPVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
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



bool kgl::SNPVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and phaseId() == cmp_snp->phaseId()
         and offset() == cmp_snp->offset()
         and variantType() == cmp_snp->variantType()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


// Order variant types.
bool kgl::SNPVariant::lessThan(const Variant& cmp_var) const {


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

  auto cmp_snp = dynamic_cast<const SNPVariant*>(&cmp_var);

  if (not cmp_snp) {

    // Must be a variant type == snp type.
    ExecEnv::log().error("SNPVariant::lessThan; Expected SNP, got: {}", cmp_var.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  if (reference() < cmp_snp->reference()) {

    return true;

  } else if (reference() > cmp_snp->reference()) {

    return false;

  } else if (mutant() < cmp_snp->mutant()) {

    return true;

  }

  return false;

}


std::string kgl::SNPVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << DNA5::convertToChar(reference()) << offsetOutput(offset(), output_index);
  ss << mutantChar() << delimiter;

  return ss.str();

}


bool kgl::SNPVariant::mutateSequence(SignedOffset_t offset_adjust,
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

  // Check the reference.
  if (reference() != dna_sequence_ptr->at(sequence_offset)) {

    ExecEnv::log().warn("mutateSequence(), SNP reference base: {} does not match sequence base: {} at genome: {} contig: {} offset: {}",
                        DNA5::convertToChar(reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        genomeId(), contigId(), offset());

  }
  // Mutate the sequence
  dna_sequence_ptr->modifyBase(sequence_offset, mutant());

  // SNPs do not modify the sequence size.
  sequence_size_modify = 0;

  return true;

}
