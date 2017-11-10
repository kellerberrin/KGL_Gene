//
// Created by kellerberrin on 31/10/17.
//

#include "kgl_variant_snp.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SNPVariant - SNPs generated from the SAM/BAM read data.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::SNPVariantDNA5::equivalent(const Variant& cmp_var) const {

  auto cmp_snp = dynamic_cast<const SNPVariantDNA5*>(&cmp_var);

  if (cmp_snp == nullptr) return false;

  return contigId() == cmp_snp->contigId()
         and contigOffset() == cmp_snp->contigOffset()
         and reference() == cmp_snp->reference()
         and mutant() == cmp_snp->mutant();

}


std::string kgl::SNPVariantDNA5::output(char delimiter, VariantOutputIndex output_index) const
{
  std::stringstream ss;
  ss << genomeOutput(delimiter, output_index);
  ss << mutation(delimiter, output_index);
  ss << mutantCount() << "/" << readCount() << delimiter;
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << NucleotideColumn_DNA5::offsetToNucleotide(idx) << ":" << countArray()[idx] << delimiter;
  }
  ss << '\n';

  return ss.str();

}


std::string kgl::SNPVariantDNA5::mutation(char delimiter, VariantOutputIndex output_index) const
{

  std::stringstream ss;

  if (not codingSequences().empty()) {

    std::shared_ptr<const CodingSequence> sequence = codingSequences().getFirst();

    ss << sequence->getGene()->id() << delimiter << sequence->getCDSParent()->id() << delimiter;

    if (NucleotideColumn_DNA5::isBaseCode(mutant())) {

      ContigOffset_t codon_offset;
      typename AminoAcidTypes::AminoType reference_amino;
      typename AminoAcidTypes::AminoType mutant_amino;

      if (contig()->SNPMutation(sequence,
                                contigOffset(),
                                reference(),
                                mutant(),
                                codon_offset,
                                reference_amino,
                                mutant_amino)) {

        ss << reference_amino << offsetOutput(codon_offset, output_index) << mutant_amino << delimiter;
        ss << reference() << offsetOutput(contigOffset(), output_index) << mutant() << delimiter;
//        sequence->getGene()->recusivelyPrintsubfeatures();



      }

    } else {  // is a deletion or insert

      ContigOffset_t sequence_offset;
      ContigSize_t sequence_length;
      DNA5Sequence::offsetWithinSequence(sequence, contigOffset(), sequence_offset, sequence_length);
      ss << reference() << offsetOutput(sequence_offset, output_index) << mutant() << delimiter;
      ss << reference() << offsetOutput(contigOffset(), output_index) << mutant() << delimiter;

    }

  } else if (not geneMembership().empty()) {

    std::shared_ptr<const GeneFeature> gene_ptr = geneMembership().front();
    ss << gene_ptr->id() << " ";
    ss << reference() << offsetOutput(contigOffset(), output_index) << mutant() << delimiter;

  } else { // non coding (non-gene) variant or unknown

    ss << reference() << offsetOutput(contigOffset(), output_index) << mutant() << delimiter;

  }

  return ss.str();

}

