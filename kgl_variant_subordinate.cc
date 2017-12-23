//
// Created by kellerberrin on 23/12/17.
//

#include "kgl_variant_subordinate.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  SubordinateSNP Variant. A virtual class held in compound variants, produces a modified text output.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::SubordinateSNP::suboutput(char delimiter, VariantOutputIndex output_index) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << quality() << delimiter;
  ss << subname() << delimiter << size() << delimiter;
  ss << submutation(delimiter, output_index);
  ss << evidence()->output(delimiter, output_index);
  ss << '\n';

  return ss.str();

}

std::string kgl::SubordinateSNP::submutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    ContigSize_t base_in_codon;
    ContigOffset_t codon_offset;

    codonOffset(codon_offset, base_in_codon);

    ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR
       << offsetOutput(base_in_codon, output_index) << delimiter;
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << ExtendDNA5::convertToChar(mutant()) << delimiter;


  } else if (not geneMembership().empty()) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << " ";
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

  } else { // non coding (non-gene) variant or unknown

    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

  }

  return ss.str();

}


// complement base if -ve strand and coding or intron.
kgl::CodingDNA5::Alphabet kgl::SubordinateSNP::strandNucleotide(DNA5::Alphabet nucleotide) const {

  switch(type()) {

    case VariantSequenceType::CDS_CODING:

      if (codingSequences().empty()) {

        ExecEnv::log().error("strandReference(), empty coding sequence for coding variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true));
        return DNA5::convertToCodingDN5(nucleotide);
      }
      switch(codingSequences().getFirst()->getCDSParent()->sequence().strand()) {

        case StrandSense::UNKNOWN:
          ExecEnv::log().error("strandReference(), unknown coding sequence for variant: {}",
                               output(' ', VariantOutputIndex::START_0_BASED, true));
          return CodingDNA5::Alphabet::N;
        case StrandSense::FORWARD:
          return DNA5::convertToCodingDN5(nucleotide);

        case StrandSense::REVERSE:
          return DNA5::complementNucleotide(nucleotide);

      }

    case VariantSequenceType::INTRON:

      if (geneMembership().empty()) {

        ExecEnv::log().error("strandReference(), no gene for intron variant: {}",
                             output(' ', VariantOutputIndex::START_0_BASED, true));
        return CodingDNA5::Alphabet::N;

      }
      switch(geneMembership().front()->sequence().strand()) {

        case StrandSense::UNKNOWN:
          ExecEnv::log().error("strandReference(), unknown coding sequence for variant: {}",
                               output(' ', VariantOutputIndex::START_0_BASED, true));
          return CodingDNA5::Alphabet::N;
        case StrandSense::FORWARD:
          return DNA5::convertToCodingDN5(nucleotide);

        case StrandSense::REVERSE:
          return DNA5::complementNucleotide(reference());

      }

    case VariantSequenceType::NON_CODING:
      return DNA5::convertToCodingDN5(nucleotide);

  }

  return DNA5::convertToCodingDN5(nucleotide); // To stop compiler complaints.

}


bool kgl::SubordinateSNP::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SubordinateSNP*>(&cmp_var);

  if (not cmp_snp) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and type() == cmp_snp->type()
         and variantType() == cmp_snp->variantType()
         and codingSequenceId() == cmp_snp->codingSequenceId()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


bool kgl::SubordinateSNP::codonMutation(ContigOffset_t &codon_offset,
                                        ContigSize_t &base_in_codon,
                                        AminoAcid::Alphabet &reference_amino,
                                        AminoAcid::Alphabet &mutant_amino) const {


  if (not codonOffset(codon_offset, base_in_codon)) {

    ExecEnv::log().error("codonMutation() called for non coding variant: {}",
                         output(' ', VariantOutputIndex::START_0_BASED, true));
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

  if (ExtendDNA5::isBaseCode(mutant())) {

    reference_amino = contig()->getAminoAcid(codon);
    codon.modifyBase(base_in_codon, strandNucleotide(ExtendDNA5::extendToBase(mutant())));
    mutant_amino = contig()->getAminoAcid(codon);


  } else {

    mutant_amino = contig()->getAminoAcid(codon);
    mutant_amino = AminoAcid::AMINO_UNKNOWN;

  }

  return true;

}



std::string kgl::SubordinateSNP::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    if (variantType() == VariantType::SNP) {

      ContigOffset_t codon_offset;
      ContigSize_t base_in_codon;
      AminoAcid::Alphabet reference_amino;
      AminoAcid::Alphabet mutant_amino;

      if (codonMutation(codon_offset, base_in_codon, reference_amino, mutant_amino)) {

        ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
        ss << offsetOutput(base_in_codon, output_index) << delimiter;
        ss << AminoAcid::convertToChar(reference_amino) << offsetOutput(codon_offset, output_index);
        ss << AminoAcid::convertToChar(mutant_amino) << delimiter;
        ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
        ss << ExtendDNA5::convertToChar(mutant()) << delimiter;
//        sequence->getGene()->recusivelyPrintsubfeatures();

      }

    } else if (variantType()== VariantType::INSERT) {  // is an insertion

      ContigSize_t base_in_codon;
      ContigOffset_t codon_offset;

      codonOffset(codon_offset, base_in_codon);

      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << "+(" << size() << ")";
      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
      ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

    } else if (variantType()== VariantType::DELETE) {  // is a deletion

      ContigSize_t base_in_codon;
      ContigOffset_t codon_offset;

      codonOffset(codon_offset, base_in_codon);

      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << "-(" << size() << ")";
      ss << offsetOutput(codon_offset, output_index) << CODON_BASE_SEPARATOR;
      ss << offsetOutput(base_in_codon, output_index) << delimiter;
      ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
      ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

    }

  } else if (type() == VariantSequenceType::INTRON) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << delimiter;
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

  } else { // else non coding (non-gene) variant or unknown

    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << ExtendDNA5::convertToChar(mutant()) << delimiter;

  }

  return ss.str();

}
