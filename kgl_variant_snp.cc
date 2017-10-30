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
  ss << genomeOutput();
  ss << reference() << (contigOffset() + 1) << mutant() << " ";
  ss << mutantCount() << "/" << readCount() << " [";
  for (size_t idx = 0; idx < countArray().size(); ++idx) {
    ss << NucleotideColumn_DNA5::offsetToNucleotide(idx) << ":" << countArray()[idx] << " ";
  }
  ss << "]";
  return ss.str();
}

