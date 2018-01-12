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
  ss << quality() << delimiter;
  ss << name() << delimiter << size() << delimiter;
  ss << mutation(delimiter, output_index);

  if (detail) {

    ss << evidence()->output(delimiter, output_index);

  }

  ss << '\n';

  return ss.str();

}



bool kgl::SNPVariant::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariant*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and type() == cmp_snp->type()
         and variantType() == cmp_snp->variantType()
         and codingSequenceId() == cmp_snp->codingSequenceId()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


bool kgl::SNPVariant::codonMutation(ContigOffset_t &codon_offset,
                                    ContigSize_t &base_in_codon,
                                    AminoAcid::Alphabet &reference_amino,
                                    AminoAcid::Alphabet &mutant_amino) const {


  if (not codonOffset(codon_offset, base_in_codon)) {

    ExecEnv::log().error("codonMutation() called for non coding variant");
    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    return false;

  }

  auto sequence_offset = static_cast<ContigOffset_t>(codon_offset * Codon::CODON_SIZE);

  const std::shared_ptr<const CodingSequence> coding_seq_ptr = codingSequences().getFirst();
  std::shared_ptr<DNA5SequenceCoding> codon_sequence = contig()->sequence().codingSubSequence(coding_seq_ptr,
                                                                                              sequence_offset,
                                                                                              Codon::CODON_SIZE);
  if (codon_sequence->length() != Codon::CODON_SIZE) {

    ExecEnv::log().error("codonMutation(), expected codon sequence size 3: got size: {}", codon_sequence->length());
    reference_amino = AminoAcid::AMINO_UNKNOWN;  // The unknown amino acid
    mutant_amino = AminoAcid::AMINO_UNKNOWN;
    return false;

  }

  Codon codon(codon_sequence, 0);  // Create the codon.

  reference_amino = contig()->getAminoAcid(codon);
  codon.modifyBase(base_in_codon, strandMutant());
  mutant_amino = contig()->getAminoAcid(codon);

  return true;

}


std::string kgl::SNPVariant::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigOffset_t codon_offset;
    ContigSize_t base_in_codon;
    AminoAcid::Alphabet reference_amino;
    AminoAcid::Alphabet mutant_amino;

    if (codonMutation(codon_offset, base_in_codon, reference_amino, mutant_amino)) {

      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << AminoAcid::convertToChar(reference_amino) << offsetOutput(codon_offset, output_index);
      ss << AminoAcid::convertToChar(mutant_amino) << delimiter;

      ContigOffset_t coding_sequence_offset;
      ContigSize_t coding_sequence_length;
      SequenceOffset::refOffsetWithinCodingSequence(sequence, offset(), coding_sequence_offset, coding_sequence_length);

      ss << CodingDNA5::convertToChar(strandReference()) << offsetOutput(coding_sequence_offset, output_index);
      ss << CodingDNA5::convertToChar(strandMutant()) << delimiter;
//        sequence->getGene()->recusivelyPrintsubfeatures();

    }


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


bool kgl::SNPVariant::mutateSequence(SignedOffset_t offset_adjust,
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

    ExecEnv::log().warn("mutateSequence(), SNP reference base: {} does not match sequence base: {} at contig: {} offset: {}",
                        DNA5::convertToChar(reference()),
                        DNA5::convertToChar(dna_sequence_ptr->at(sequence_offset)),
                        contig()->contigId(), offset());

  }
  // Mutate the sequence
  dna_sequence_ptr->modifyBase(sequence_offset, mutant());

  return true;

}
