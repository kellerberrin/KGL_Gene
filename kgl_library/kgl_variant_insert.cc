//
// Created by kellerberrin on 23/12/17.
//


#include "kgl_variant_single.h"
#include "kgl_sequence_offset.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DeleteVariant - generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::InsertVariant::output(char delimiter, VariantOutputIndex output_index, bool detail) const
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


bool kgl::InsertVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const InsertVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and phaseId() == cmp_snp->phaseId()
         and offset() == cmp_snp->offset()
         and variantType() == cmp_snp->variantType()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


// Order variant types.
bool kgl::InsertVariant::lessThan(const Variant& cmp_var) const {


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

  auto cmp_insert = dynamic_cast<const InsertVariant*>(&cmp_var);

  if (not cmp_insert) {

    // Must be a variant type == insert type.
    ExecEnv::log().error("InsertVariant::lessThan; Expected InsertVariant, got: {}", cmp_var.output(' ', VariantOutputIndex::START_0_BASED, false));
    return false;

  }

  if (reference() < cmp_insert->reference()) {

    return true;

  } else if (reference() > cmp_insert->reference()) {

    return false;

  } else if (mutant() < cmp_insert->mutant()) {

    return true;

  }

  return false;

}


std::string kgl::InsertVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  ss << DNA5::convertToChar(reference()) << offsetOutput(offset(), output_index);
  ss << mutantChar() << delimiter;

  return ss.str();

}



bool kgl::InsertVariant::mutateSequence(SignedOffset_t offset_adjust,
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


    ExecEnv::log().info("mutateSequence(), Insert reference base: {} does not match sequence base: {}; Genome: {} Contig: {} Offset: {}",
                        DNA5::convertToChar(reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        genomeId(), contigId(), offset());

  }
  // Mutate the sequence
  StringDNA5 seq_string;
  seq_string.push_back(mutant());
  DNA5SequenceLinear insert_seq(seq_string);

  if (dna_sequence_ptr->insertSubSequence(sequence_offset, insert_seq)) {

    sequence_size_modify = insert_seq.length();

  } else {

    sequence_size_modify = 0;

  }

  return true;

}
