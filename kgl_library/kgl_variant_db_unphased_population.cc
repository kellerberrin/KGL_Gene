//
// Created by kellerberrin on 13/08/18.
//


#include "kgl_variant_db_unphased_population.h"

#include <fstream>
#include <vector>
#include <string>


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

  // Create an array of genomes
  for (const auto& [genome, genome_ptr] : getMap()) {

      auto genome_copy = genome_ptr->deepCopy();

      if (not population_copy->addGenome(genome_copy)) {

        ExecEnv::log().error( "UnphasedPopulation::deepCopy, could not insert genome: {} (duplicate) into population: {}",
                              genome, population_copy->populationId());

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
    auto insert_result = genome_map_.try_emplace(genome_id, genome_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedPopulation::getCreateGenome(), Could not add genome: {} to the population: {}", genome_id, populationId());
      return std::nullopt;

    }

    return genome_ptr;

  }

}


bool kgl::UnphasedPopulation::addGenome(const std::shared_ptr<UnphasedGenome>& genome_ptr) {

  auto result = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  if (not result.second) {

    ExecEnv::log().error("UnphasedPopulation::addGenome(), could not add genome: {} (duplicate) to the population", genome_ptr->genomeId());

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

      ExecEnv::log().critical("UnphasedPopulation::filterVariants(), could not add filtered genome: {} to the population", genome_id);

    }

  }

  return filtered_population_ptr;

}


bool kgl::UnphasedPopulation::addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                          const std::vector<GenomeId_t>& genome_vector) {

  bool result = true;
  for (auto& genome : genome_vector) {

    std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add/create genome: {}", genome);
      result = false;
      continue;

    }

    if (not genome_opt.value()->addVariant(variant_ptr)) {

      ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add variant to genome: {}", genome);
      result = false;
      continue;

    }

  }

  return result;

}

// Only adds the variant if it does not already exist.
bool kgl::UnphasedPopulation::addUniqueVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                                const std::vector<GenomeId_t>& genome_vector) {

  bool result = true;
  for (auto& genome : genome_vector) {

    std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add/create genome: {}", genome);
      result = false;
      continue;

    }

    bool result = genome_opt.value()->addUniqueVariant(variant_ptr);
    if (not result) {

      ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add variant to genome: {}", genome);
      result = false;
      continue;

    }

  }

  return result;

}


// Unconditional Merge.
size_t kgl::UnphasedPopulation::mergePopulation(const std::shared_ptr<const UnphasedPopulation>& merge_population) {

  size_t variant_count = 0;
  bool result = true;
  for (const auto& [genome, genome_ptr] : merge_population->getMap()) {

    std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add/create genome: {}", genome);
      result = false;
      continue;
    }

    variant_count += genome_opt.value()->mergeGenome(genome_ptr);

  }

  return variant_count;

}


std::pair<size_t, size_t> kgl::UnphasedPopulation::validate(const std::shared_ptr<const GenomeReference>& genome_db) const {

  std::pair<size_t, size_t> population_count{0, 0};
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::pair<size_t, size_t> genome_count = genome_ptr->validate(genome_db);

    if (genome_count.first != genome_count.second) {

      ExecEnv::log().warn("UnphasedPopulation::validate(), Population: {} Failed to Validate Genome: {}, Total Variants: {}, Validated: {}",
                    populationId(), genome_id, genome_count.first, genome_count.second);

    }

    population_count.first += genome_count.first;
    population_count.second += genome_count.second;

  }

  return population_count;

}



std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedPopulation::compressPopulation() const {

  std::shared_ptr<kgl::UnphasedGenome> compressedGenome(std::make_shared<UnphasedGenome>("Compressed"));

  if (not processAll(*compressedGenome, &UnphasedGenome::addVariant)) {

    ExecEnv::log().error("UnphasedPopulation::compressPopulation(); problem compressing population: {}", populationId());
  }

  return  compressedGenome;

}


std::optional<std::shared_ptr<const kgl::InfoEvidenceHeader>> kgl::UnphasedPopulation::getVCFInfoEvidenceHeader() const {

  for (auto const& genome : getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second) {

          if (variant_ptr->evidence().infoData()) {

            return variant_ptr->evidence().infoData().value()->evidenceHeader();

          }

        }

      }

    }

  }

  return std::nullopt;

}


std::shared_ptr<kgl::UnphasedPopulation> kgl::UnphasedPopulation::uniqueVariantPopulation() const {

  // Create a unique variant (no duplicate) version of this population.
  std::shared_ptr<UnphasedPopulation>unique_ptr(std::make_shared<UnphasedPopulation>(populationId()));

  for (auto const& [genome_id, genome_ptr] : getMap()) {

    auto unique_genome = genome_ptr->uniqueGenome();

    if (not unique_ptr->addGenome(unique_genome)) {

      ExecEnv::log().error("UnphasedPopulation::uniqueVariantPopulation, problem adding unique genome: {}", genome_id) ;

    }

  }

  return unique_ptr;

}



