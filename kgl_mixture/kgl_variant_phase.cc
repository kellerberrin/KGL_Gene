//
// Created by kellerberrin on 23/04/18.
//

#include "kgl_variant_phase.h"
#include "kgl_variant_evidence.h"
#include "kgl_variant_single.h"
#include <fstream>


namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Not Thread safe.
// This object accepts unphased variants and trivally phases them as haploid.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


// Just copy into a population object.
bool kgl::GenomePhasing::haploidPhasing(std::shared_ptr<const UnphasedPopulation> unphased_population_ptr,
                                        std::shared_ptr<const GenomeDatabase> genome_db,
                                        std::shared_ptr<PhasedPopulation> haploid_population)  {

  bool result = true;

  for (auto genome : unphased_population_ptr->getMap()) {

    // Create the GenomeVariant object.
    std::shared_ptr<GenomeVariant> genome_variant = GenomeVariant::emptyGenomeVariant(genome.first, GenomeVariant::HAPLOID_GENOME, genome_db);

    // Add all the contig.
    for (auto contig : genome.second->getMap()) {

      // Add all offsets.
      for (auto offset : contig.second->getMap()) {

        // add all the variants
        for (auto variant : offset.second) {

          std::shared_ptr<Variant> mutable_variant = std::const_pointer_cast<Variant>(variant);
          mutable_variant->updatePhaseId(ContigVariant::HAPLOID_HOMOLOGOUS_INDEX);   // Assign to the first (and only) homologous contig.
          genome_variant->addVariant(variant);

        } // All variants.

      } // All offsets

    } // All contigs

    if (not haploid_population->addGenomeVariant(genome_variant)) {

      ExecEnv::log().error("Unable to add genome: {} to haploid population", genome_variant->genomeId());
      result = false;

    }

  } // All genomes

  ExecEnv::log().info("Haploid Phasing, pre_phase variants: {}, genomes: {}, resultant population variants: {}, genomes: {}",
                      unphased_population_ptr->variantCount(), unphased_population_ptr->getMap().size(),
                      haploid_population->variantCount(), haploid_population->getMap().size());

  return result;

}



