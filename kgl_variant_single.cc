//
// Created by kellerberrin on 23/12/17.
//


#include "kgl_variant_single.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  SingleVariant Variant. A virtual class held in compound variants, produces a modified text output.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::string kgl::SingleVariant::suboutput(char delimiter, VariantOutputIndex output_index) const
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

std::string kgl::SingleVariant::submutation(char delimiter, VariantOutputIndex output_index) const
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
    ss << mutantChar() << delimiter;


  } else if (not geneMembership().empty()) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << " ";
    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  } else { // non coding (non-gene) variant or unknown

    ss << DNA5::convertToChar(reference()) << offsetOutput(contigOffset(), output_index);
    ss << mutantChar() << delimiter;

  }

  return ss.str();

}


// complement base if -ve strand and coding or intron.
kgl::CodingDNA5::Alphabet kgl::SingleVariant::strandNucleotide(DNA5::Alphabet nucleotide) const {

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



