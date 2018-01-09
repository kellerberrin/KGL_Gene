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
  ss << quality() << delimiter;
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence()->output(delimiter, output_index);

  }

  ss << '\n';

  return ss.str();

}


bool kgl::InsertVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const InsertVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and type() == cmp_snp->type()
         and variantType() == cmp_snp->variantType()
         and codingSequenceId() == cmp_snp->codingSequenceId()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


std::string kgl::InsertVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigSize_t base_in_codon;
    ContigOffset_t codon_offset;

    codonOffset(codon_offset, base_in_codon);

    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
    ss << offsetOutput(base_in_codon, output_index) << delimiter;
    ss << "+(" << size() << ")";
    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
    ss << offsetOutput(base_in_codon, output_index) << delimiter;

    ContigOffset_t coding_sequence_offset;
    ContigSize_t coding_sequence_length;
    SequenceOffset::refOffsetWithinCodingSequence(sequence, offset(), coding_sequence_offset, coding_sequence_length);

    ss << CodingDNA5::convertToChar(strandReference()) << offsetOutput(coding_sequence_offset, output_index);
    ss << CodingDNA5::convertToChar(strandMutant()) << delimiter;

  } else if (type() == VariantSequenceType::INTRON) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << delimiter;
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  } else { // else non coding (non-gene) variant or unknown

    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  }

  return ss.str();

}



bool kgl::InsertVariant::mutateSequence(SignedOffset_t offset_adjust,
                                        std::shared_ptr<DNA5SequenceLinear> dna_sequence_ptr) const {


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

    ExecEnv::log().info("insertSequence(), Insert reference base: {} does not match sequence base: {} at insert contig offset: {}",
                        DNA5::convertToChar(reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        offset());

  }
  // Mutate the sequence
  StringDNA5 seq_string;
  seq_string.push_back(mutant());
  DNA5SequenceLinear insert_seq(seq_string);
  dna_sequence_ptr->insertSubSequence(sequence_offset, insert_seq);

  return true;

}
