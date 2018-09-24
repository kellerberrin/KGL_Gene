//
// Created by kellerberrin on 23/04/18.
//


#include <memory>
#include <fstream>
#include "kgl_patterns.h"
#include "kgl_variant_db.h"
#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container to hold genome variants for populations
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::PhasedPopulation::getCreateGenome(const GenomeId_t& genome_id,
                                            PhaseId_t ploidy,
                                            const std::shared_ptr<const GenomeDatabase>& genome_db,
                                            std::shared_ptr<GenomeVariant>& genome) {

  auto result = population_variant_map_.find(genome_id);

  if (result != population_variant_map_.end()) {

    genome = std::const_pointer_cast<GenomeVariant>(result->second);

    if (genome->ploidy() != ploidy) {

      ExecEnv::log().error("PhasedPopulation::getCreateGenome(), Genome: {} Ploidy: {} does not equal Requested Ploidy: {}",
                           genome_id, genome->ploidy(), ploidy);
      return false;

    }

    return true;

  } else {

    genome = GenomeVariant::emptyGenomeVariant(genome_id, ploidy, genome_db);

    std::pair<GenomeId_t, std::shared_ptr<const GenomeVariant>> insert_genome(genome->genomeId(), genome);

    auto result = population_variant_map_.insert(insert_genome);

    if (not result.second) {

      ExecEnv::log().critical("PhasedPopulation::getCreateGenome(), Serious Error, could not add genome: {} to the population", genome_id);

    }

    return result.second;

  }

}



bool kgl::PhasedPopulation::getGenomeVariant(const GenomeId_t& genome_id,
                                             std::shared_ptr<const GenomeVariant>& genome_variant) const {

  auto result = population_variant_map_.find(genome_id);

  if (result != population_variant_map_.end()) {

    genome_variant = result->second;
    return true;

  } else {

    genome_variant = nullptr;
    return false;

  }

}


bool kgl::PhasedPopulation::addGenomeVariant(std::shared_ptr<const GenomeVariant> genome_variant) {

  auto result = population_variant_map_.insert(std::pair<GenomeId_t, std::shared_ptr<const GenomeVariant>>(genome_variant->genomeId(), genome_variant));

  return result.second;

}


size_t kgl::PhasedPopulation::variantCount() const {

  size_t variant_count = 0;
  for (auto genome : getMap()) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


std::shared_ptr<kgl::PhasedPopulation> kgl::PhasedPopulation::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::PhasedPopulation> filtered_population_ptr(std::make_shared<kgl::PhasedPopulation>(populationId()));

  for (const auto& genome_variant : population_variant_map_) {

    std::shared_ptr<kgl::GenomeVariant> filtered_genome_ptr = genome_variant.second->filterVariants(filter);
    filtered_population_ptr->addGenomeVariant(filtered_genome_ptr);
    ExecEnv::log().vinfo("Genome: {} has: {} filtered variants", genome_variant.first, filtered_genome_ptr->variantCount());

  }

  return filtered_population_ptr;

}


std::shared_ptr<kgl::PhasedPopulation> kgl::PhasedPopulation::filterGenomes(const PopulationId_t& population_id,
                                                                            const std::vector<GenomeId_t>& list) const {

  std::shared_ptr<kgl::PhasedPopulation> filtered_population_ptr(std::make_shared<kgl::PhasedPopulation>(population_id));

  for (const auto& genome_variant : population_variant_map_) {

    for (auto genome : list) {

      if (genome == genome_variant.first) {

        if (not filtered_population_ptr->addGenomeVariant(genome_variant.second->deepCopy())) {

          ExecEnv::log().error("Filtered Population: {} unable to add genome: {} (duplicate)",
                              filtered_population_ptr->populationId(), genome_variant.first);

        }

        break;

      }

    }

  }

  ExecEnv::log().info("Filtered Population: {} has: {} filtered genomes with: {} variants",
                      filtered_population_ptr->populationId(), filtered_population_ptr->getMap().size(), filtered_population_ptr->variantCount());

  return filtered_population_ptr;

}
