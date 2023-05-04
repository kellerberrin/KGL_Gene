//
// Created by kellerberrin on 4/05/23.
//


#include "kgl_variant_db_population.h"
#include "kgl_variant_filter.h"

namespace kgl = kellerberrin::genome;




// Multi-tasked filtering for large populations.
std::unique_ptr<kgl::PopulationDB> kgl::PopulationDB::viewFilter(const BaseFilter& filter) const {

  // Create the new population.
  std::unique_ptr<PopulationDB> filtered_population_ptr(std::make_unique<PopulationDB>(populationId(), dataSource()));

  // Edge Condition, if no genomes then simply exit.
  if (getMap().empty()) {

    return filtered_population_ptr;

  }

  // Only a population filter is implemented at this level.
  std::shared_ptr<const FilterPopulations> population_filter = std::dynamic_pointer_cast<const FilterPopulations>(filter.clone());
  if (population_filter) {

    return population_filter->applyFilter(*this);

  }

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), WorkflowThreads::defaultThreads());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<std::shared_ptr<GenomeDB>>> future_vector;
  // Required by the thread pool.
  auto filter_lambda = [](std::shared_ptr<const GenomeDB> genome_ptr,
                          std::shared_ptr<const BaseFilter> filter_ptr)-> std::shared_ptr<GenomeDB> {

    return genome_ptr->viewFilter(*filter_ptr);

  };

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::shared_ptr<const BaseFilter> filter_ptr = filter.clone();
    std::future<std::shared_ptr<GenomeDB>> future = thread_pool.enqueueFuture(filter_lambda, genome_ptr, filter_ptr);
    future_vector.push_back(std::move(future));

  }


  // Add in the filtered genomes.
  for (auto& future : future_vector) {

    std::shared_ptr<GenomeDB> filtered_genome_ptr = future.get();
    if (not filtered_population_ptr->addGenome(filtered_genome_ptr)) {

      ExecEnv::log().error("PopulationDB::filter; could not add filtered genome to the population");

    }

  }

  return filtered_population_ptr;

}

// Multi-threaded filtering for large populations.
// We can do this because smart pointer reference counting (only) is thread safe.
// Returns a std::pair with .first the original number of variants, .second the filtered number of variants.
std::pair<size_t, size_t> kgl::PopulationDB::selfFilter(const BaseFilter& filter) {

  // This routine modifies the populationDB data structure, so only permit one thread at a time.
  // Note the routine is internally multi-threaded.
  std::scoped_lock lock(insitufilter_mutex_);

  // Edge Condition, if no genomes then simply exit.
  if (getMap().empty()) {

    return {0, 0};

  }

  // Only a population filter is implemented at this level.
  std::shared_ptr<const FilterPopulations> population_filter = std::dynamic_pointer_cast<const FilterPopulations>(filter.clone());
  if (population_filter) {

    size_t prior_count = variantCount();

    auto population_ptr = population_filter->applyFilter(*this);
    genome_map_ = std::move(population_ptr->genome_map_);

    size_t post_count = variantCount();

    return {prior_count, post_count};

  }

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), WorkflowThreads::defaultThreads());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;
  // Required by the thread pool.

  auto filter_lambda = [](std::shared_ptr<GenomeDB> genome_ptr,
                          std::shared_ptr<const BaseFilter> filter_ptr) -> std::pair<size_t, size_t> {

    return genome_ptr->selfFilter(*filter_ptr);

  };

  // Queue a thread for each genome.
  for (auto& [genome_id, genome_ptr] : getMap()) {

    std::shared_ptr<const BaseFilter> filter_ptr = filter.clone();
    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueFuture(filter_lambda, genome_ptr, filter_ptr);
    future_vector.push_back(std::move(future));

  }

  // Wait for the threads to finish.
  std::pair<size_t, size_t> filter_counts{0, 0};
  for (auto& future : future_vector) {

    auto [original_count, filtered_count] = future.get();
    filter_counts.first += original_count;
    filter_counts.second += filtered_count;

  }

  return filter_counts;

}


