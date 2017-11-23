//
// Created by kellerberrin on 22/11/17.
//



#include "kgl_variant_factory.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The compound variant factory.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<kgl::GenomeVariant>
kgl::VariantCompoundFactory::disaggregateCompoundVariants(const std::shared_ptr<const GenomeVariant>& genome_variant,
                                                          const std::shared_ptr<const GenomeDatabase>& genome_db) const {

  std::shared_ptr<kgl::GenomeVariant> disaggreagated = genome_variant->emptyGenomeVariant(genome_variant->genomeId(),
                                                                                          genome_db);

  for (auto contig_variant : genome_variant->contigMap()) {

    for (auto variant : contig_variant.second->getMap()) {

      if (variant.second->isCompound()) {

        std::shared_ptr<const CompoundVariant>
        compound_variant = std::static_pointer_cast<const CompoundVariant>(variant.second);

        for (auto single_variant : compound_variant->getMap()) {

          if (not disaggreagated->addVariant(single_variant.second)) {

            ExecEnv::log().error("Cannot add disaggregated variant: {} - same contig offset as existing variant",
                                 single_variant.second->output(' ', VariantOutputIndex::START_0_BASED));

          }

        }

      }

    }

  }

  return disaggreagated;

}


