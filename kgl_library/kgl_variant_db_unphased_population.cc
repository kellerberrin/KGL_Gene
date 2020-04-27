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


std::shared_ptr<kgl::UnphasedPopulation> kgl::UnphasedPopulation::deepCopy() const {

  // Duplicate population ids are not a problem.
  std::shared_ptr<UnphasedPopulation> population_copy(std::make_shared<UnphasedPopulation>(populationId()));

  for (auto const& [genome_id, genome_ptr] : getMap()) {

    if (not population_copy->addGenome(genome_ptr->deepCopy())) {

      ExecEnv::log().critical("UnphasedPopulation::deepCopy(), probable duplicate, could not add genome: {} to the population", genome_id);

    }

  }

  return population_copy;

}


std::optional<std::shared_ptr<kgl::UnphasedGenome>> kgl::UnphasedPopulation::getCreateGenome(const GenomeId_t& genome_id) {

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    std::shared_ptr<UnphasedGenome> genome_ptr = std::make_shared<UnphasedGenome>(genome_id);
    std::pair<GenomeId_t, std::shared_ptr<UnphasedGenome>> new_genome(genome_id, genome_ptr);
    auto result = genome_map_.insert(new_genome);

    if (not result.second) {

      ExecEnv::log().error("UnphasedPopulation::getCreateGenome(), Could not add genome: {} to the population: {}", genome_id, populationId());
      return std::nullopt;

    }

    return genome_ptr;

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

  for (auto const& [contig_id, contig_ptr] : result->second->getMap()) {

    for (auto const& [offset, variant_vector] : contig_ptr->getMap()) {

      if (variant_vector.empty()) {

        ExecEnv::log().error(
        "UnphasedPopulation::genomePhasingStats(); Zero sized variant vector, genome: {}, contig: {}, offset: {}",
        genome_id, contig_id, offset);
        return_result = false;
        continue;

      }

      switch(variant_vector.size()) {

        case 0: {
          ExecEnv::log().error(
          "UnphasedPopulation::genomePhasingStats(); Zero sized variant vector, genome: {}, contig: {}, offset: {}",
          genome_id, contig_id, offset);
          return_result = false;
        }
          break;

        case 1:  ++singleheterozygous;
          break;

        case 2: {

          if (variant_vector[0]->equivalent(*variant_vector[1])) {

            ++homozygous;

          } else {

            ++heterozygous;

          }

        }
          break;

        default:
          ExecEnv::log().warn("UnphasedPopulation::genomePhasingStats(); {} variants found at genome: {}, contig: {}, offset: {}",
                              variant_vector.size(), genome_id, contig_id, offset);


      } // switch

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

  std::shared_ptr<kgl::UnphasedPopulation> filtered_population_ptr(std::make_shared<kgl::UnphasedPopulation>(populationId()));

  for (const auto& [genome_id, genome_ptr] : getMap()) {

    std::shared_ptr<kgl::UnphasedGenome> filtered_genome_ptr = genome_ptr->filterVariants(filter);
    if (not filtered_population_ptr->addGenome(filtered_genome_ptr)) {

      ExecEnv::log().critical("UnphasedPopulation::filterVariants(), could not add filtered genome: {} to the population", genome_ptr->genomeId());

    }

    ExecEnv::log().vinfo("Genome: {} has: {} filtered variants", genome_id, filtered_genome_ptr->variantCount());

  }

  return filtered_population_ptr;

}


bool kgl::UnphasedPopulation::addVariant(std::shared_ptr<const Variant>& variant_ptr) {

  std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(variant_ptr->genomeId());
  if (not genome_opt) {

    ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add/create genome: {}", variant_ptr->genomeId());
    return false;

  }

  if (not genome_opt.value()->addVariant(variant_ptr)) { // thread safe

    ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add variant to genome: {}", variant_ptr->genomeId());
    return false;

  }

  return true;

}


void kgl::UnphasedPopulation::mergePopulation(std::shared_ptr<const UnphasedPopulation> merge_population) {

  for (auto const& genome : merge_population->getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& [offset, variant_vector] : contig.second->getMap()) {

        for (auto variant_ptr : variant_vector) {

          if (not addVariant(variant_ptr)) {

            ExecEnv::log().error("UnphasedPopulation::mergePopulation(); Cannot merge variant offset: {} from: {} to to: {}",
                                  offset, merge_population->populationId(), populationId());

          }

        }

      }

    }

  }

}


bool kgl::UnphasedPopulation::validate(const std::shared_ptr<const GenomeDatabase>& genome_db) const {

  bool result = true;
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    if (not genome_ptr->validate(genome_db)) {

      result = false;
      ExecEnv::log().warn("UnphasedPopulation::validate(), Population: {} failed to validate Variants in Genome: {}", populationId(), genome_id);

    }

  }

  return result;

}


