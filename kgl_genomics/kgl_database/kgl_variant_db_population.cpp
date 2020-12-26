//
// Created by kellerberrin on 25/12/20.
//

#include "kgl_variant_db_population.h"

namespace kgl = kellerberrin::genome;


// This function is used by VCF parsers to create a variant database.
// The function has been made thread safe for multiple parser thread access.
std::optional<std::shared_ptr<kgl::GenomeVariantArray>> kgl::PopulationVariant::getCreateGenome(const GenomeId_t& genome_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    std::shared_ptr<GenomeVariantArray> genome_ptr = std::make_shared<GenomeVariantArray>(genome_id);
    auto insert_result = genome_map_.try_emplace(genome_id, genome_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("UnphasedPopulation::getCreateGenome(), Could not add genome: {} to the population: {}", genome_id, populationId());
      return std::nullopt;

    }

    return genome_ptr;

  }

}


std::optional<std::shared_ptr<kgl::GenomeVariantArray>> kgl::PopulationVariant::getGenome(const GenomeId_t& genome_id) const {


  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


bool kgl::PopulationVariant::addGenome(const std::shared_ptr<GenomeVariantArray>& genome_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  if (not result.second) {

    ExecEnv::log().error("UnphasedPopulation::addGenome(), could not add genome: {} (duplicate) to the population", genome_ptr->genomeId());

  }

  return result.second;

}


size_t kgl::PopulationVariant::variantCount() const {

  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), ThreadPool::hardwareThreads());
  ThreadPool thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<size_t>> future_vector;
  // Required by the thread pool.
  /// todo: This could be re-coded as a lambda, investigate threadpool type deduction for lambda functions.
  struct CountClass {

    static size_t countGenome(std::shared_ptr<const GenomeVariantArray> genome_ptr) { return genome_ptr->variantCount(); };

  } ;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::future<size_t> future = thread_pool.enqueueTask(&CountClass::countGenome, genome_ptr);
    future_vector.push_back(std::move(future));

  }


  size_t variant_count = 0;

  // Add the genome variant counts.
  for (auto& future : future_vector) {

    variant_count += future.get();

  }

  return variant_count;

}




// Multi-tasking filtering for large populations.
// We can do this because smart pointer reference counting (only) is thread safe.
std::shared_ptr<kgl::PopulationVariant> kgl::PopulationVariant::filterVariants(const VariantFilter& filter) const {

  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), ThreadPool::hardwareThreads());
  ThreadPool thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<std::shared_ptr<GenomeVariantArray>>> future_vector;
  // Required by the thread pool.
  /// todo: This could be re-coded as a lambda, investigate threadpool type deduction for lambda functions.
  std::shared_ptr<const VariantFilter> filter_ptr = filter.clone();
  struct FilterClass {

    static std::shared_ptr<GenomeVariantArray>
    filterGenome(std::shared_ptr<GenomeVariantArray> genome_ptr,
                 std::shared_ptr<const VariantFilter>& filter_ptr) { return genome_ptr->filterVariants(*filter_ptr); };


  } ;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::future<std::shared_ptr<GenomeVariantArray>> future = thread_pool.enqueueTask(&FilterClass::filterGenome, genome_ptr, filter_ptr);
    future_vector.push_back(std::move(future));

  }

  // Create the new population.
  std::shared_ptr<PopulationVariant> filtered_population_ptr(std::make_shared<PopulationVariant>(populationId(), dataSource()));

  // Add in the filtered genomes.
  for (auto& future : future_vector) {

    std::shared_ptr<GenomeVariantArray> filtered_genome_ptr = future.get();
    if (not filtered_population_ptr->addGenome(filtered_genome_ptr)) {

      ExecEnv::log().error("PopulationVariant::filterVariants; could not add filtered genome to the population");

    }

  }

  return filtered_population_ptr;

}

// Multi-threaded filtering for large populations.
// We can do this because smart pointer reference counting (only) is thread safe.
// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::PopulationVariant::inSituFilter(const VariantFilter& filter) {

  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), ThreadPool::hardwareThreads());
  ThreadPool thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;
  // Required by the thread pool.
  /// todo: This could be re-coded as a lambda, investigate threadpool type deduction for lambda functions.
  std::shared_ptr<const VariantFilter> filter_ptr = filter.clone();
  struct FilterClass {

    static std::pair<size_t, size_t>
    inSituFilterGenome(std::shared_ptr<GenomeVariantArray> genome_ptr,
                       std::shared_ptr<const VariantFilter>& filter_ptr) { return genome_ptr->inSituFilter(*filter_ptr); };

  } ;

  // Queue a thread for each genome.
  for (auto& [genome_id, genome_ptr] : getMap()) {

    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueTask(&FilterClass::inSituFilterGenome, genome_ptr, filter_ptr);
    future_vector.push_back(std::move(future));

  }

  // Wait for the threads to finish.
  std::pair<size_t, size_t> filter_counts{0, 0};
  for (auto& future : future_vector) {

    auto genome_count = future.get();
    filter_counts.first += genome_count.first;
    filter_counts.second += genome_count.second;

  }

  // Delete empty genomes.
  for (auto it = genome_map_.begin(); it != genome_map_.end(); ++it) {

    if (it->second->getMap().empty()) {

      it = genome_map_.erase(it);

    }

  }

  return filter_counts;

}




bool kgl::PopulationVariant::addVariant( const std::shared_ptr<const Variant>& variant_ptr,
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
size_t kgl::PopulationVariant::mergePopulation(const std::shared_ptr<const PopulationVariant>& merge_population) {

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
std::pair<size_t, size_t> kgl::PopulationVariant::validate(const std::shared_ptr<const GenomeReference>& genome_db) const {

  ThreadPool thread_pool(ThreadPool::hardwareThreads());
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    // function, object_ptr, arg1
    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueTask(&GenomeVariantArray::validate,
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


std::shared_ptr<kgl::GenomeVariantArray> kgl::PopulationVariant::compressPopulation() const {

  std::shared_ptr<GenomeVariantArray> compressedGenome(std::make_shared<GenomeVariantArray>("Compressed"));

  if (not processAll(*compressedGenome, &GenomeVariantArray::addVariant)) {

    ExecEnv::log().error("UnphasedPopulation::compressPopulation(); problem compressing population: {}", populationId());
  }

  return  compressedGenome;

}


std::shared_ptr<kgl::GenomeVariantArray> kgl::PopulationVariant::uniqueUnphasedGenome() const {

  std::shared_ptr<GenomeVariantArray> unphasedGenome(std::make_shared<GenomeVariantArray>("UniqueUnphased"));

  if (not processAll(*unphasedGenome, &GenomeVariantArray::addUniqueUnphasedVariant)) {

    ExecEnv::log().error("UnphasedPopulation::UniqueUnphased(); problem creating unique unphased genome with population: {}", populationId());
  }

  return  unphasedGenome;

}


std::optional<std::shared_ptr<const kgl::InfoEvidenceHeader>> kgl::PopulationVariant::getVCFInfoEvidenceHeader() const {

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

