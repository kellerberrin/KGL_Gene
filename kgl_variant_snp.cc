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


std::string kgl::SNPVariantDNA5::output() const
{
  std::stringstream ss;
  ss << genomeOutput() << " ";
  ss << mutation();
  ss << mutantCount() << "/" << readCount() << " [";
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << NucleotideColumn_DNA5::offsetToNucleotide(idx) << ":" << countArray()[idx] << " ";
  }
  ss << "]" << "\n";

  return ss.str();

}


std::string kgl::SNPVariantDNA5::mutation() const
{

  std::stringstream ss;

  if (type() == VariantSequenceType::CDS_CODING) {

    std::shared_ptr<const CodingSequence> sequence;
    codingSequences(sequence);

    ss << sequence->getGene()->id() << " " << sequence->getCDSParent()->id() << " ";

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

        ss << reference_amino << codon_offset << mutant_amino << " ";
//        sequence->getGene()->recusivelyPrintsubfeatures();

      ss << reference() << contigOffset() << mutant() << " ";


      } else {  // is a deletion or insert

        ContigOffset_t sequence_offset;
        ContigSize_t sequence_length;
        DNA5Sequence::offsetWithinSequence(sequence, contigOffset(), sequence_offset, sequence_length);
        ss << reference() << sequence_offset << mutant() << " ";
        ss << reference() << contigOffset() << mutant() << " ";

      }


    }

  } else if (type() == VariantSequenceType::INTRON) {

    std::shared_ptr<const GeneFeature> gene_ptr;
   if (not geneMembership(gene_ptr)) {

     ExecEnv::log().error("Intron does not have a valid gene");

   } else {

     ss << gene_ptr->id() << " ";

   }

    ss << reference() << contigOffset() << mutant() << " ";


  } else { // non coding (non-gene) variant or unknown

    ss << reference() << contigOffset() << mutant() << " ";

  }

  return ss.str();

}

