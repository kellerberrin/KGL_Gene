//
// Created by kellerberrin on 23/04/23.
//

#include "kgl_variant_filter_db_genome.h"
#include "kgl_variant_filter_db_contig.h"

namespace kgl = kellerberrin::genome;


/// Filters a population to the listed genomes (if they exist).
/// Note that this is a shallow copy of the reference population.
/// Use selfFilter() or deepCopy() to create a permanent population view.
std::unique_ptr<kgl::PopulationDB> kgl::GenomeListFilter::applyFilter(const PopulationDB& population) const {

  std::unique_ptr<PopulationDB> filtered_population_ptr(std::make_unique<PopulationDB>(population.populationId(), population.dataSource()));

  for (auto const& genome_id : genome_set_) {

    auto iter = population.getMap().find(genome_id);
    if (iter != population.getMap().end()) {

      if (not filtered_population_ptr->addGenome(iter->second)) {

        ExecEnv::log().warn("GenomeListFilter::applyFilter; Unexpected duplicate genome: {} in population: {}", genome_id, population.populationId());

      }

    }

  }

  return filtered_population_ptr;

}
