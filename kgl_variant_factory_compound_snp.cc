//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory_compound.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound SNP variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Compound SNPs are created to hold 2 or 3 coding SNPs that modify an Amino Acid simultaneously.
// 1. The SNP can only be in the same coding sequence (same mRNA parent).
// 2. They must be within an offset of 2 of each other.
// 3. They must lie within the same codon boundary (a subset of point 2).



std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantCompoundSNPFactory::compoundSNP(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                            const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<GenomeVariant>
  genome_compoundSNP_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_variants->genomeId(), genome_db_ptr);


  return genome_compoundSNP_variants;

}


