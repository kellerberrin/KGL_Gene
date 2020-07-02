//
// Created by kellerberrin on 13/08/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_POPULATION_H
#define KGL_VARIANT_DB_UNPHASED_POPULATION_H


#include "kgl_variant_db_genome.h"

#include <map>
#include <mutex>
#include <functional>

namespace kellerberrin::genome {   //  organization::project


////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// An internal parser variant object that holds variants until they can be phased.
// This object hold variants for a population.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class VariantGenome>
using VariantGenomeMap = std::map<ContigId_t, std::shared_ptr<VariantGenome>>;

template<class VariantGenome>
class PopulationVariant;

using UnphasedPopulation = PopulationVariant<UnphasedGenome>;

template<class VariantGenome>
class PopulationVariant {

public:

  explicit PopulationVariant(const PopulationId_t& population_id) : population_id_(population_id) {}
  PopulationVariant(const PopulationVariant&) = delete; // Use deep copy.
  virtual ~PopulationVariant() = default;

  PopulationVariant& operator=(const PopulationVariant&) = delete; // Use deep copy.

  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<PopulationVariant> deepCopy() const;

  // Create the genome variant if it does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<VariantGenome>> getCreateGenome(const GenomeId_t& genome_id);

  // Retrieve a genome
  [[nodiscard]] std::optional<std::shared_ptr<VariantGenome>> getGenome(const GenomeId_t& genome_id) const;

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] std::vector<GenomeId_t> genomeList() const;

  [[nodiscard]] std::shared_ptr<PopulationVariant> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const VariantGenomeMap<VariantGenome>& getMap() const { return genome_map_; }

  void clear() { genome_map_.clear(); }

  // Generate phasing statistics (only valid with ploidy 2 variants).
  [[nodiscard]] bool genomePhasingStats( const GenomeId_t& genome_id,
                                         size_t& heterozygous,
                                         size_t& homozygous,
                                         size_t& singleheterozygous) const;

  [[nodiscard]] bool addGenome(const std::shared_ptr<VariantGenome>& genome);

  // Unconditionally add a variant to the population.
  [[nodiscard]] bool addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                 const std::vector<GenomeId_t>& genome_vector);

  // The first bool is normal operation. The second bool is if a unique variant was added to the population.
  [[nodiscard]] bool addUniqueVariant( const std::shared_ptr<const Variant>& variant,
                                       const std::vector<GenomeId_t>& genome_vector);

  [[nodiscard]] const PopulationId_t& populationId() const { return population_id_; }
  void setPopulationId(const PopulationId_t& population_id) { population_id_ = population_id; }
  // unconditionally merge (retains duplicates) genomes and variants into this population.
  [[nodiscard]] size_t mergePopulation(const std::shared_ptr<const PopulationVariant>& merge_population);
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the reference genome.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db) const;
  // Compress a population into a single genome. Done when generating aggregate variant statistics for a population.
  [[nodiscard]] std::shared_ptr<VariantGenome> compressPopulation() const;
  // Get the Info header, get the field header object from the first variant in the population.
  // Careful, this implicitly assumes that all variants in the population have the same DataMemoryBlock this only true
  // Of populations generated by a single VCF file.
  [[nodiscard]] std::optional<std::shared_ptr<const InfoEvidenceHeader>> getVCFInfoEvidenceHeader() const;
  // Processes all variants in the population with class Obj and Func = &Obj::objFunc(const shared_ptr<const Variant>&)
  template<class Obj, typename Func> bool processAll(Obj& object, Func objFunc) const;
  // Create a population of unique variants. All duplicate variants are removed.
  [[nodiscard]] std::shared_ptr<PopulationVariant> uniqueVariantPopulation() const;


private:

  VariantGenomeMap<VariantGenome> genome_map_;
  PopulationId_t population_id_;
  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

};

// General purpose population processing template.
// Processes all variants in the population with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class VariantGenome>
template<class Obj, typename Func>
bool PopulationVariant<VariantGenome>::processAll(Obj& object, Func objFunc)  const {

  for (auto const& [genome, genome_ptr] : getMap()) {

    if (not genome_ptr->processAll(object, objFunc)) {

      ExecEnv::log().error("UnphasedPopulation::processAll<Obj, Func>; error with genome: {}", genome);
      return false;

    }

  }

  return true;

}


// This function should always be used to copy a variant database.
template<class VariantGenome>
std::shared_ptr<PopulationVariant<VariantGenome>> PopulationVariant<VariantGenome>::deepCopy() const {

  // Duplicate population ids are not a problem.
  std::shared_ptr<PopulationVariant> population_copy(std::make_shared<PopulationVariant>(populationId()));

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

// This function is used by VCF parsers to create a variant database.
// The function has been made thread safe for multiple parser thread access.
template<class VariantGenome>
std::optional<std::shared_ptr<VariantGenome>> PopulationVariant<VariantGenome>::getCreateGenome(const GenomeId_t& genome_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    std::shared_ptr<VariantGenome> genome_ptr = std::make_shared<VariantGenome>(genome_id);
    auto insert_result = genome_map_.try_emplace(genome_id, genome_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedPopulation::getCreateGenome(), Could not add genome: {} to the population: {}", genome_id, populationId());
      return std::nullopt;

    }

    return genome_ptr;

  }

}


template<class VariantGenome>
std::optional<std::shared_ptr<VariantGenome>> PopulationVariant<VariantGenome>::getGenome(const GenomeId_t& genome_id) const {


  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


template<class VariantGenome>
bool PopulationVariant<VariantGenome>::addGenome(const std::shared_ptr<VariantGenome>& genome_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  if (not result.second) {

    ExecEnv::log().error("UnphasedPopulation::addGenome(), could not add genome: {} (duplicate) to the population", genome_ptr->genomeId());

  }

  return result.second;

}


template<class VariantGenome>
size_t PopulationVariant<VariantGenome>::variantCount() const {

  size_t variant_count = 0;

  for (auto genome : genome_map_) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


template<class VariantGenome>
bool PopulationVariant<VariantGenome>::genomePhasingStats( const GenomeId_t& genome_id,
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

      if (variant_vector->getVariantArray().empty()) {

        ExecEnv::log().error(
        "UnphasedPopulation::genomePhasingStats(); Zero sized variant vector, genome: {}, contig: {}, offset: {}",
        genome_id, contig_id, offset);
        return_result = false;
        continue;

      }

      switch(variant_vector->getVariantArray().size()) {

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

          if (variant_vector->getVariantArray()[0]->equivalent(*variant_vector->getVariantArray()[1])) {

            ++homozygous;

          } else {

            ++heterozygous;

          }

        }
          break;

        default:
          ExecEnv::log().warn("UnphasedPopulation::genomePhasingStats(); {} variants found at genome: {}, contig: {}, offset: {}",
                              variant_vector->getVariantArray().size(), genome_id, contig_id, offset);


      } // switch

    } // offset

  } // contig

  return return_result;

}

template<class VariantGenome>
std::vector<GenomeId_t> PopulationVariant<VariantGenome>::genomeList() const {

  std::vector<GenomeId_t> genome_list;

  for (auto genome : getMap()) {

    genome_list.push_back(genome.first);

  }

  return genome_list;

}

template<class VariantGenome>
std::shared_ptr<PopulationVariant<VariantGenome>> PopulationVariant<VariantGenome>::filterVariants(const VariantFilter& filter) const {

  std::shared_ptr<PopulationVariant> filtered_population_ptr(std::make_shared<PopulationVariant>(populationId()));

  for (const auto& [genome_id, genome_ptr] : getMap()) {

    std::shared_ptr<VariantGenome> filtered_genome_ptr = genome_ptr->filterVariants(filter);
    if (not filtered_population_ptr->addGenome(filtered_genome_ptr)) {

      ExecEnv::log().error("PopulationVariant::filterVariants(), could not add filtered genome: {} (duplicate) to the population", genome_id);

    }

  }

  return filtered_population_ptr;

}

template<class VariantGenome>
bool PopulationVariant<VariantGenome>::addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                                   const std::vector<GenomeId_t>& genome_vector) {

  bool result = true;
  for (auto& genome : genome_vector) {

    auto genome_opt = getCreateGenome(genome);
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


// Unconditional Merge.
template<class VariantGenome>
size_t PopulationVariant<VariantGenome>::mergePopulation(const std::shared_ptr<const PopulationVariant>& merge_population) {

  size_t variant_count = 0;
  for (const auto& [genome, genome_ptr] : merge_population->getMap()) {

    auto genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("UnphasedPopulation::addUniqueVariant; Could not add/create genome: {}", genome);
      continue;
    }

    variant_count += genome_opt.value()->mergeGenome(genome_ptr);

  }

  return variant_count;

}


template<class VariantGenome>
std::pair<size_t, size_t> PopulationVariant<VariantGenome>::validate(const std::shared_ptr<const GenomeReference>& genome_db) const {

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


template<class VariantGenome>
std::shared_ptr<VariantGenome> PopulationVariant<VariantGenome>::compressPopulation() const {

  std::shared_ptr<VariantGenome> compressedGenome(std::make_shared<VariantGenome>("Compressed"));

  if (not processAll(*compressedGenome, &VariantGenome::addVariant)) {

    ExecEnv::log().error("UnphasedPopulation::compressPopulation(); problem compressing population: {}", populationId());
  }

  return  compressedGenome;

}


template<class VariantGenome>
std::optional<std::shared_ptr<const InfoEvidenceHeader>> PopulationVariant<VariantGenome>::getVCFInfoEvidenceHeader() const {

  for (auto const& genome : getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second->getVariantArray()) {

          if (variant_ptr->evidence().infoData()) {

            return variant_ptr->evidence().infoData().value()->evidenceHeader();

          }

        }

      }

    }

  }

  return std::nullopt;

}

template<class VariantGenome>
std::shared_ptr<PopulationVariant<VariantGenome>> PopulationVariant<VariantGenome>::uniqueVariantPopulation() const {

  // Create a unique variant (no duplicate) version of this population.
  std::shared_ptr<PopulationVariant> unique_pop_ptr(std::make_shared<PopulationVariant>(populationId()));

  for (auto const& [genome_id, genome_ptr] : getMap()) {

    auto unique_genome = genome_ptr->uniqueGenome();

    if (not unique_pop_ptr->addGenome(unique_genome)) {

      ExecEnv::log().error("UnphasedPopulation::uniqueVariantPopulation, problem adding unique genome: {}", genome_id) ;

    }

  }

  return unique_pop_ptr;

}



}   // end namespace




#endif //KGL_VARIANT_DB_UNPHASED_POPULATION_H
