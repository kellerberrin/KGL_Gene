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

      ExecEnv::log().critical("UnphasedPopulation::filterVariants(), could not add filtered genome: {} to the population", genome_id);

    }

  }

  return filtered_population_ptr;

}



// Test if an equivalent variant already exists in the genome.
bool kgl::UnphasedPopulation::variantExists(std::shared_ptr<const Variant> variant) const {

  auto result = getMap().find (variant->contigId());
  if (result == getMap().end()) {

    return false;

  }

  return result->second->variantExists(variant);

}


bool kgl::UnphasedPopulation::addVariant(std::shared_ptr<const Variant> variant_ptr) {

  std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(variant_ptr->genomeId());
  if (not genome_opt) {

    ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add/create genome: {}", variant_ptr->genomeId());
    return false;

  }

  if (not genome_opt.value()->addVariant(variant_ptr)) {

    ExecEnv::log().error("UnphasedPopulation::addVariant; Could not add variant to genome: {}", variant_ptr->genomeId());
    return false;

  }

  return true;

}

// Only adds the variant if it does not already exist.
bool kgl::UnphasedPopulation::addUniqueVariant(std::shared_ptr<const Variant> variant_ptr) {

  std::optional<std::shared_ptr<UnphasedGenome>> genome_opt = getCreateGenome(variant_ptr->genomeId());
  if (not genome_opt) {

    ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add/create genome: {}", variant_ptr->genomeId());
    return false;

  }

  bool result = genome_opt.value()->addUniqueVariant(variant_ptr);
  if (not result) {

    ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add variant to genome: {}", variant_ptr->genomeId());
    return false;

  }

  return true;

}


// Unconditional Merge.
void kgl::UnphasedPopulation::mergePopulation(std::shared_ptr<const UnphasedPopulation> merge_population) {

  if (merge_population->processAll(*this, &UnphasedPopulation::addVariant)) {

    ExecEnv::log().error("UnphasedPopulation::mergePopulation(); Cannot merge population: {} with population: {}",
                          merge_population->populationId(), populationId());

  }

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


std::pair<size_t, size_t> kgl::UnphasedPopulation::mergeUniqueGenome(const std::shared_ptr<const UnphasedGenome> genome) {

  std::pair<size_t, size_t> merge_count{0, 0};

  for (auto const& contig : genome->getMap()) {

    for (auto const& [offset, variant_vector] : contig.second->getMap()) {

      for (auto const& variant_ptr : variant_vector) {

        ++merge_count.first;
        bool result = addUniqueVariant(variant_ptr);
        if (not result) {

          ExecEnv::log().error("UnphasedPopulation::mergeUniqueGenome(); Cannot merge variant offset: {} from Genome: {} to Population: {}",
                               offset, genome->genomeId(), populationId());

        }

      }

    }

  }

  merge_count.second = variantCount();

  return merge_count;

}


std::pair<size_t, size_t> kgl::UnphasedPopulation::mergeUniquePopulation(const std::shared_ptr<const UnphasedPopulation> merged_population) {

  std::pair<size_t, size_t> merge_count{0, 0};

  for (auto const& genome : merged_population->getMap()) {

    std::pair<size_t, size_t> result = mergeUniqueGenome(genome.second);
    merge_count.first += result.first;
    merge_count.second += result.second;

  }

  return merge_count;

}

std::shared_ptr<kgl::UnphasedGenome> kgl::UnphasedPopulation::compressPopulation() const {

  std::shared_ptr<kgl::UnphasedGenome> compressedGenome(std::make_shared<UnphasedGenome>("Compressed"));

  for (auto const& genome : getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second) {

          if (not compressedGenome->addVariant(variant_ptr)) {

            ExecEnv::log().error("UnphasedPopulation::compressPopulation(); Cannot Add variant: {}",
            variant_ptr->output(' ', VariantOutputIndex::START_0_BASED, false));

          }

        }

      }

    }

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



std::shared_ptr<const kgl::UnphasedPopulation> kgl::UnphasedPopulation::setUnion(std::shared_ptr<const UnphasedPopulation> pop_b) const {

  // Create a no duplicate variants version of this population.
  std::shared_ptr<UnphasedPopulation>union_ptr(std::make_shared<UnphasedPopulation>("Union(" + populationId() + ", " + pop_b->populationId() + ")"));

  // All all unique variants to the union population.
  if (not processAll(*union_ptr, &UnphasedPopulation::addUniqueVariant)) {

    ExecEnv::log().error("UnphasedPopulation::setUnion, problem adding unique variants from population: {}", populationId());

  }

  // Add all unique variants from the other ('B') population.
  if (not pop_b->processAll(*union_ptr, &UnphasedPopulation::addUniqueVariant)) {

    ExecEnv::log().error("UnphasedPopulation::setUnion, problem adding unique variants from population: {}", pop_b->populationId());

  }

  return union_ptr;

}


std::shared_ptr<const kgl::UnphasedPopulation> kgl::UnphasedPopulation::setIntersection(std::shared_ptr<const UnphasedPopulation> pop_b) const {

  std::shared_ptr<UnphasedPopulation>intersect_ptr(std::make_shared<UnphasedPopulation>("Intersection(" + populationId() + ", " + pop_b->populationId() + ")"));

  // Temp Local class to perform the operation.
  struct Intersection {

    std::shared_ptr<UnphasedPopulation> intersect_ptr_;
    std::shared_ptr<const UnphasedPopulation> population_b_ptr_;

    bool addIntersection(std::shared_ptr<const Variant> variant_ptr) {

      if (population_b_ptr_->variantExists(variant_ptr)) {

        return intersect_ptr_->addUniqueVariant(variant_ptr);

      }

      return true;

    }

  };

  Intersection intersect;
  intersect.intersect_ptr_ = intersect_ptr;
  intersect.population_b_ptr_ = pop_b;

  // Add all unique variants common to 'A' and 'B' populations.
  if (not processAll(intersect, &Intersection::addIntersection)) {

    ExecEnv::log().error("UnphasedPopulation::setIntersection, problem adding unique variants from population: {}", pop_b->populationId());

  }

  return intersect_ptr;

}


std::shared_ptr<const kgl::UnphasedPopulation> kgl::UnphasedPopulation::setComplement(std::shared_ptr<const UnphasedPopulation> pop_b) const {

  std::shared_ptr<UnphasedPopulation>complement_ptr(std::make_shared<UnphasedPopulation>(populationId() + " complement(" + pop_b->populationId() + ")"));

  // Temp Local class to perform the operation.
  struct Complement {

    std::shared_ptr<UnphasedPopulation> complement_ptr_;
    std::shared_ptr<const UnphasedPopulation> population_b_ptr_;

    bool addComplement(std::shared_ptr<const Variant> variant_ptr) {

      if (not population_b_ptr_->variantExists(variant_ptr)) {

        return complement_ptr_->addUniqueVariant(variant_ptr);

      }

      return true;

    }

  };

  Complement complement;
  complement.complement_ptr_ = complement_ptr;
  complement.population_b_ptr_ = pop_b;

  // Add all unique variants in 'A' but not in 'B' populations.
  if (not processAll(complement, &Complement::addComplement)) {

    ExecEnv::log().error("UnphasedPopulation::setComplement, problem complementing unique variants from population: {}", pop_b->populationId());

  }

  return complement_ptr;

}



