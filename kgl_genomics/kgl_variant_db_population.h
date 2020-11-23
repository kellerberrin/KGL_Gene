//
// Created by kellerberrin on 13/08/18.
//

#ifndef KGL_VARIANT_DB_UNPHASED_POPULATION_H
#define KGL_VARIANT_DB_UNPHASED_POPULATION_H


#include "kgl_variant_db_genome.h"
#include "kgl_data_base.h"
#include "kel_thread_pool.h"


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

template<class VariantGenome, class PopulationBase>
class PopulationVariant : public PopulationBase {

public:

  explicit PopulationVariant(const PopulationId_t& population_id) : PopulationBase(population_id) {}
  PopulationVariant(const PopulationVariant&) = delete; // Use deep copy.
  ~PopulationVariant() override = default;

  PopulationVariant& operator=(const PopulationVariant&) = delete; // Use deep copy.

  [[nodiscard]] const PopulationId_t& populationId() const { return PopulationBase::Id(); }

  void setPopulationId(const PopulationId_t& population_id) { PopulationBase::setId(population_id); }
  // Use this to copy the object.
  [[nodiscard]] std::shared_ptr<PopulationVariant> deepCopy() const;

  // Create the genome variant if it does not exist.
  [[nodiscard]] std::optional<std::shared_ptr<VariantGenome>> getCreateGenome(const GenomeId_t& genome_id);

  // Retrieve a genome
  [[nodiscard]] std::optional<std::shared_ptr<VariantGenome>> getGenome(const GenomeId_t& genome_id) const;

  [[nodiscard]] size_t variantCount() const;

  [[nodiscard]] std::shared_ptr<PopulationVariant> filterVariants(const VariantFilter& filter) const;

  [[nodiscard]] const VariantGenomeMap<VariantGenome>& getMap() const { return genome_map_; }

  void clear() { genome_map_.clear(); }

  bool addGenome(const std::shared_ptr<VariantGenome>& genome);

  // Unconditionally add a variant to the population.
  // This function is thread safe for concurrent updates.
  // The population structure cannot be 'read' while it is being updated.
  [[nodiscard]] bool addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                 const std::vector<GenomeId_t>& genome_vector);

  // unconditionally merge (retains duplicates) genomes and variants into this population.
  [[nodiscard]] size_t mergePopulation(const std::shared_ptr<const PopulationVariant>& merge_population);
  // Validate returns a pair<size_t, size_t>. The first integer is the number of variants examined.
  // The second integer is the number variants that pass inspection by comparison to the reference genome.
  [[nodiscard]] std::pair<size_t, size_t> validate(const std::shared_ptr<const GenomeReference>& genome_db) const;
  // Compress a population into a single genome. Done when generating aggregate variant statistics for a population.
  [[nodiscard]] std::shared_ptr<VariantGenome> compressPopulation() const;
  // Compress a population into a single genome of unique (only) variants. Removes any variant phasing information.
  // Used to convert a Diploid/Haploid population to an unphased single genome. Source populations are unchanged.
  template<class UnphasedGenome> [[nodiscard]] std::shared_ptr<UnphasedGenome> uniqueUnphasedGenome() const;
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
  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

};

// General purpose population processing template.
// Processes all variants in the population with class Obj and Func = &(bool Obj::objFunc(const std::shared_ptr<const Variant>))
template<class VariantGenome, class PopulationBase>
template<class Obj, typename Func>
bool PopulationVariant<VariantGenome, PopulationBase>::processAll(Obj& object, Func objFunc)  const {

  for (auto const& [genome, genome_ptr] : getMap()) {

    if (not genome_ptr->processAll(object, objFunc)) {

      ExecEnv::log().error("UnphasedPopulation::processAll<Obj, Func>; error with genome: {}", genome);
      return false;

    }

  }

  return true;

}


// This function should always be used to copy a variant database.
template<class VariantGenome, class PopulationBase>
std::shared_ptr<PopulationVariant<VariantGenome, PopulationBase>> PopulationVariant<VariantGenome, PopulationBase>::deepCopy() const {

  // Duplicate population ids are not a problem.
  std::shared_ptr<PopulationVariant> population_copy(std::make_shared<PopulationVariant>(PopulationBase::populationId()));

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
template<class VariantGenome, class PopulationBase>
std::optional<std::shared_ptr<VariantGenome>> PopulationVariant<VariantGenome, PopulationBase>::getCreateGenome(const GenomeId_t& genome_id) {

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


template<class VariantGenome, class PopulationBase>
std::optional<std::shared_ptr<VariantGenome>> PopulationVariant<VariantGenome, PopulationBase>::getGenome(const GenomeId_t& genome_id) const {


  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


template<class VariantGenome, class PopulationBase>
bool PopulationVariant<VariantGenome, PopulationBase>::addGenome(const std::shared_ptr<VariantGenome>& genome_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  if (not result.second) {

    ExecEnv::log().error("UnphasedPopulation::addGenome(), could not add genome: {} (duplicate) to the population", genome_ptr->genomeId());

  }

  return result.second;

}


template<class VariantGenome, class PopulationBase>
size_t PopulationVariant<VariantGenome, PopulationBase>::variantCount() const {

  size_t variant_count = 0;

  for (auto const& genome : genome_map_) {

    variant_count += genome.second->variantCount();

  }

  return variant_count;

}


template<class VariantGenome, class PopulationBase>
std::shared_ptr<PopulationVariant<VariantGenome, PopulationBase>> PopulationVariant<VariantGenome, PopulationBase>::filterVariants(const VariantFilter& filter) const {

  std::shared_ptr<PopulationVariant> filtered_population_ptr(std::make_shared<PopulationVariant>(populationId()));

  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::shared_ptr<VariantGenome> filtered_genome_ptr = genome_ptr->filterVariants(filter);
    if (not filtered_population_ptr->addGenome(filtered_genome_ptr)) {

      ExecEnv::log().error("PopulationVariant::filterVariants(), could not add filtered genome: {} (duplicate) to the population", genome_id);

    }

  }

  return filtered_population_ptr;

}

template<class VariantGenome, class PopulationBase>
bool PopulationVariant<VariantGenome, PopulationBase>::addVariant( const std::shared_ptr<const Variant>& variant_ptr,
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
template<class VariantGenome, class PopulationBase>
size_t PopulationVariant<VariantGenome, PopulationBase>::mergePopulation(const std::shared_ptr<const PopulationVariant>& merge_population) {

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


// Ensures that all variants are correctly specified.
template<class VariantGenome, class PopulationBase>
std::pair<size_t, size_t> PopulationVariant<VariantGenome, PopulationBase>::validate(const std::shared_ptr<const GenomeReference>& genome_db) const {

  ThreadPool thread_pool;
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    // function, object_ptr, arg1
    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueTask(&VariantGenome::validate,
                                                                            genome_ptr,
                                                                            genome_db);
    future_vector.push_back(std::move(future));

  }

  // Check the results of the validation.
  std::pair<size_t, size_t> population_count{0, 0};
  for (auto& future : future_vector) {

    std::pair<size_t, size_t> genome_count = future.get();

    if (genome_count.first != genome_count.second) {

      ExecEnv::log().warn("UnphasedPopulation::validate(), Population: {} Failed to Validate, Total Variants: {}, Validated: {}",
                          populationId(), genome_count.first, genome_count.second);

    }

    population_count.first += genome_count.first;
    population_count.second += genome_count.second;

  }

  return population_count;

}


template<class VariantGenome, class PopulationBase>
std::shared_ptr<VariantGenome> PopulationVariant<VariantGenome, PopulationBase>::compressPopulation() const {

  std::shared_ptr<VariantGenome> compressedGenome(std::make_shared<VariantGenome>("Compressed"));

  if (not processAll(*compressedGenome, &VariantGenome::addVariant)) {

    ExecEnv::log().error("UnphasedPopulation::compressPopulation(); problem compressing population: {}", populationId());
  }

  return  compressedGenome;

}



template<class VariantGenome, class PopulationBase>
template<class UnphasedGenome>
std::shared_ptr<UnphasedGenome> PopulationVariant<VariantGenome, PopulationBase>::uniqueUnphasedGenome() const {

  std::shared_ptr<UnphasedGenome> unphasedGenome(std::make_shared<UnphasedGenome>("UniqueUnphased"));

  if (not processAll(*unphasedGenome, &UnphasedGenome::addUniqueUnphasedVariant)) {

    ExecEnv::log().error("UnphasedPopulation::UniqueUnphased(); problem creating unique unphased genome with population: {}", populationId());
  }

  return  unphasedGenome;

}


template<class VariantGenome, class PopulationBase>
std::optional<std::shared_ptr<const InfoEvidenceHeader>> PopulationVariant<VariantGenome, PopulationBase>::getVCFInfoEvidenceHeader() const {

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


}   // end namespace




#endif //KGL_VARIANT_DB_UNPHASED_POPULATION_H
