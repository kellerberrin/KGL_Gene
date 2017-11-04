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
  ss << mutation() << " ";
  ss << mutantCount() << "/" << readCount() << " [";
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << NucleotideColumn_DNA5::offsetToNucleotide(idx) << ":" << countArray()[idx] << " ";
  }
  ss << "]";
  return ss.str();
}


std::string kgl::SNPVariantDNA5::mutation() const
{

  std::stringstream ss;
  if (type() != VariantSequenceType::NON_CODING) {

    ss << "Gene(s):";
    GeneVector gene_vector = geneMembership();

    for (auto gene : gene_vector) {

      ss << gene->id() << " ";

      if (type() == VariantSequenceType::CDS_CODING) {

        SortedCDSVector sorted_cds_vec;
        gene->getSortedCDS(sorted_cds_vec);

        for (auto sorted_cds : sorted_cds_vec) {

          if (NucleotideColumn_DNA5::isBaseCode(mutant())) {

            ContigOffset_t codon_offset;
            typename AminoAcidTypes::AminoType reference_amino;
            typename AminoAcidTypes::AminoType mutant_amino;

            contig()->SNPMutation(sorted_cds,
                                  contigOffset(),
                                  reference(),
                                  mutant(),
                                  codon_offset,
                                  reference_amino,
                                  mutant_amino);

            ss << reference_amino << codon_offset << mutant_amino << " ";
            ss << reference() << contigOffset() << mutant() << " ";
//          std::shared_ptr<AminoSequence> amino_sequence = contig()->getAminoSequence(sorted_cds);
//          ss << "\n" << amino_sequence->getProteinString() << "\n";
//          gene->recusivelyPrintsubfeatures(0);

          } else {

            ContigOffset_t sequence_offset;
            ContigSize_t sequence_length;
            DNA5Sequence::offsetWithinSequence(sorted_cds, contigOffset(), sequence_offset, sequence_length);
            ss << reference() << sequence_offset << mutant() << " ";
            ss << reference() << contigOffset() << mutant() << " ";

          } // if mutant is base code

        } // for all gene CDS sequences

      } else { // the variant must be an intron

        ss << reference() << contigOffset() << mutant() << " ";

      }

    } // for all genes.

  } else { // non coding (non-gene) variant

    ss << reference() << contigOffset() << mutant() << " ";

  }

  return ss.str();

}

