//
// Created by kellerberrin on 13/08/18.
//


#include "kgl_variant_db_unphased_population.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An object that holds variants until they can be phased.
// This object hold variants for a population.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::UnphasedPopulation::getCreateGenome(const GenomeId_t& genome_id,
                                              std::shared_ptr<UnphasedGenome>& genome) {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    genome = result->second;
    return true;

  } else {

    genome = std::make_shared<UnphasedGenome>(genome_id);
    std::pair<GenomeId_t, std::shared_ptr<UnphasedGenome>> new_genome(genome_id, genome);
    auto result = genome_map_.insert(new_genome);

    if (not result.second) {

      ExecEnv::log().critical("UnphasedPopulation::getCreateGenome(), Serious Error, could not add genome: {} to the population", genome_id);

    }

    return result.second;

  }

}


bool kgl::UnphasedPopulation::addGenome(std::shared_ptr<UnphasedGenome> genome_ptr) {

  std::pair<GenomeId_t, std::shared_ptr<UnphasedGenome>> add_genome(genome_ptr->genomeId(), genome_ptr);
  auto result = genome_map_.insert(add_genome);

  if (not result.second) {

    ExecEnv::log().error("UnphasedPopulation::addGenome(), could not add genome: {} to the population", genome_ptr->genomeId());

  }

  return result.second;

}



size_t kgl::UnphasedPopulation::variantCount() const {

  size_t variant_count = 0;

  for (auto genome : genome_map_) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


bool kgl::UnphasedPopulation::genomePhasingStats(const GenomeId_t& genome_id,
                                                 size_t& heterozygous,
                                                 size_t& homozygous,
                                                 size_t& singleheterozygous) const {

  bool return_result = true;

  heterozygous = 0;
  homozygous = 0;
  singleheterozygous = 0;

  auto result = genome_map_.find(genome_id);

  if (result == genome_map_.end()) {

    ExecEnv::log().error("UnphasedPopulation::genomePhasingStats(); Could not find genome: {}", genome_id);
    return false;

  }

  for (auto contig : result->second->getMap()) {

    for (auto offset : contig.second->getMap()) {

      if (offset.second.size() == 0) {

        ExecEnv::log().error(
        "UnphasedPopulation::genomePhasingStats(); Zero sized variant vector, genome: {}, contig: {}, offset: {}",
        genome_id, contig.first, offset.first);
        return_result = false;
        continue;

      }

      if (offset.second.size() == 1) {

        if (offset.second.front().second == 1) {

          ++singleheterozygous;

        } else if (offset.second.front().second == 2) {

          ++homozygous;

        } else {

          ++homozygous;
          ExecEnv::log().warn("UnphasedPopulation::genomePhasingStats(); {} copies of variant: {}",
                              offset.second.front().second,
                              offset.second.front().first->output(' ', VariantOutputIndex::START_0_BASED, false));

        }


      } else if (offset.second.size() == 2) {

        ++heterozygous;

      } else {

        ++heterozygous;
        ExecEnv::log().warn("UnphasedPopulation::genomePhasingStats(); {} variants found at genome: {}, contig: {}, offset: {}",
                            offset.second.size(), genome_id, contig.first, offset.first);

        for (auto variant : offset.second) {

          ExecEnv::log().warn("UnphasedPopulation::genomePhasingStats(); count: {}, variant {}",
                              variant.second, variant.first->output(' ', VariantOutputIndex::START_0_BASED, false));

        }

      }

    } // offset

  } // contig

  return return_result;

}



void kgl::UnphasedPopulation::popStatistics() const {

  size_t total_variants = 0;

  for (auto genome : getMap()) {

    size_t genome_variant_count = genome.second->variantCount();

    ExecEnv::log().info("Genome: {}, Unphased Variant Count:{}", genome.first, genome_variant_count);

    total_variants += genome_variant_count;

  }

  ExecEnv::log().info("Total Unphased Variant Count:{}", total_variants);

}

std::vector<kgl::GenomeId_t> kgl::UnphasedPopulation::genomeList() const {

  std::vector<kgl::GenomeId_t> genome_list;

  for (auto genome : getMap()) {

    genome_list.push_back(genome.first);

  }

  return genome_list;

}


std::shared_ptr<kgl::UnphasedPopulation> kgl::UnphasedPopulation::filterVariants(const kgl::VariantFilter& filter) const {

  std::shared_ptr<kgl::UnphasedPopulation> filtered_population_ptr(std::make_shared<kgl::UnphasedPopulation>());

  for (const auto& genome_variant : getMap()) {

    std::shared_ptr<kgl::UnphasedGenome> filtered_genome_ptr = genome_variant.second->filterVariants(filter);
    filtered_population_ptr->addGenome(filtered_genome_ptr);
    ExecEnv::log().vinfo("Genome: {} has: {} filtered variants", genome_variant.first, filtered_genome_ptr->variantCount());

  }

  return filtered_population_ptr;

}

