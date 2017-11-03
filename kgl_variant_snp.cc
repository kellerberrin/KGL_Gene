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
  if (type() == VariantSequenceType::CDS_CODING) {

    GeneVector gene_vector = geneMembership();
    SortedCDSVector sorted_cds_vec;
    ContigOffset_t sequence_offset;
    ContigSize_t sequence_size;

    for (auto gene : gene_vector) {

      gene->getSortedCDS(sorted_cds_vec);

      for (auto sorted_cds : sorted_cds_vec) {

        DNA5Sequence::offsetWithinSequence(sorted_cds, contigOffset(), sequence_offset, sequence_size);
        if (NucleotideColumn_DNA5::isBaseCode(mutant())) {

          ss << reference() << static_cast<long>(sequence_offset/3) << mutant() << " ";

        } else {

          ss << reference() << sequence_offset << mutant() << " ";

        }

      }

    }

  } else {

    ss << reference() << contigOffset() << mutant() << " ";

  }

  return ss.str();

}

