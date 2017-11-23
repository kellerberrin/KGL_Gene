//
// Created by kellerberrin on 22/11/17.
//


#include "kgl_variant_factory.h"
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
kgl::VariantCompoundSNPFactory::compoundSNP(const std::shared_ptr<const GenomeVariant>& SNPs,
                                            const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<GenomeVariant>
  genome_compoundSNP_variants = kgl::GenomeVariant::emptyGenomeVariant(SNPs->genomeId(), genome_db_ptr);

  // Filter to only coding base mutation SNPs.
  std::shared_ptr<GenomeVariant> mutant_coding_SNPs = SNPs->filterVariants(MutantSNPFilter());
  mutant_coding_SNPs = mutant_coding_SNPs->filterVariants(InCDSFilter());

  CompoundVariantMap variant_map;
  for (auto contig_variant : mutant_coding_SNPs->contigMap()) {

    for (auto variant : contig_variant.second->getMap()) {
/*
      int64_t offset = variant.second->contigOffset();
      int64_t previous_offset = variant_map.rbegin()->second->contigOffset();
      int64_t offset_diff = ::llabs(offset - previous_offset);
      switch(variant_map.size()) {

        case 2:
          if (offset_diff == 1) {


          }

          break;

        case 1:
          break;

        case 0:
          break;

      } // switch map size
*/
    } // for variant

  } // for contig.

  return genome_compoundSNP_variants;

}


