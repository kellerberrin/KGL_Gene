//
// Created by kellerberrin on 23/11/17.
//


#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound insert variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::VariantInsertFactory::aggregateVariant(const std::shared_ptr<const Variant>& variant_ptr) const {

  if (variant_ptr->isSNP()) {

    const std::shared_ptr<const SNPVariantDNA5> SNP_ptr = std::static_pointer_cast<const SNPVariantDNA5>(variant_ptr);
    return ExtendDNA5::isInsertion(SNP_ptr->mutant());

  } else {

    return false;

  }

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantInsertFactory::compoundInsert(const std::shared_ptr<const GenomeVariant>& genome_variants,
                                          const std::shared_ptr<const GenomeDatabase>& genome_db_ptr) {

  std::shared_ptr<GenomeVariant>
  genome_compound_variants = kgl::GenomeVariant::emptyGenomeVariant(genome_variants->genomeId(), genome_db_ptr);

  return genome_compound_variants;

}
