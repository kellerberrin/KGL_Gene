//
// Created by kellerberrin on 25/12/20.
//

#include "kgl_variant_db_population.h"
#include "kgl_variant_filter.h"

namespace kgl = kellerberrin::genome;


// Use this to copy the object. Just the trivial 'TrueFilter'.
std::shared_ptr<kgl::PopulationDB> kgl::PopulationDB::deepCopy() const {

  // Can use the shallow filter because all variants are copied across
  return viewFilter(TrueFilter());

}

// Use this to empty the object. Just the trivial 'FalseFilter'.
std::pair<size_t, size_t> kgl::PopulationDB::clear() {

  auto variant_count = selfFilter(FalseFilter());
  trimEmpty();
  return variant_count;

}

// This function is used by VCF parsers to create a variant database.
// The function has been made thread safe for multiple parser thread access.
std::optional<std::shared_ptr<kgl::GenomeDB>> kgl::PopulationDB::getCreateGenome(const GenomeId_t& genome_id) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    std::shared_ptr<GenomeDB> genome_ptr = std::make_shared<GenomeDB>(genome_id);
    auto insert_result = genome_map_.try_emplace(genome_id, genome_ptr);

    if (not insert_result.second) {

      ExecEnv::log().error("PopulationDB::getCreateGenome(), Could not add genome: {} to the population: {}", genome_id, populationId());
      return std::nullopt;

    }

    return genome_ptr;

  }

}


std::optional<std::shared_ptr<kgl::GenomeDB>> kgl::PopulationDB::getGenome(const GenomeId_t& genome_id) const {


  auto result = genome_map_.find(genome_id);

  if (result != genome_map_.end()) {

    return result->second;

  } else {

    return std::nullopt;

  }

}


bool kgl::PopulationDB::addGenome(const std::shared_ptr<GenomeDB>& genome_ptr) {

  // Lock this function to concurrent access.
  std::scoped_lock lock(add_variant_mutex_);

  auto [iterator, result] = genome_map_.try_emplace(genome_ptr->genomeId(), genome_ptr);

  if (not result) {

    ExecEnv::log().error("PopulationDB::addGenome(), could not add genome: {} (duplicate) to the population", genome_ptr->genomeId());

  }

  return result;

}

// Multi-thread for speed.
size_t kgl::PopulationDB::variantCount() const {

  // Check edge condition.
  if (getMap().empty()) {

    return 0;

  }
  // Calc how many threads required.
  size_t thread_count = std::min(getMap().size(), WorkflowThreads::defaultThreads());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<size_t>> future_vector;
  // Thread pool work lambda
  auto count_lambda =  [](const std::shared_ptr<const GenomeDB>& genome_ptr)->size_t { return genome_ptr->variantCount(); };
  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    std::future<size_t> future = thread_pool.enqueueFuture(count_lambda, genome_ptr);
    future_vector.push_back(std::move(future));

  }

  size_t variant_count{0};
  // Add the genome variant counts.
  for (auto& future : future_vector) {

    variant_count += future.get();

  }

  return variant_count;

}



size_t kgl::PopulationDB::trimEmpty() {

  size_t delete_count{0};
  // Delete empty genomes.
  auto it = genome_map_.begin();
  while (it != genome_map_.end()) {

    if (it->second->variantCount() == 0) {

      it = genome_map_.erase(it);
      ++delete_count;

    } else {

      ++it;

    }

  }

  return delete_count;

}

size_t kgl::PopulationDB::squareContigs() {


  // Get a set of all contigs in all genomes.
  std::set<std::string> contig_set;
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

      contig_set.insert(contig_id);

    }

  }

  // Ensure each genome has all contigs in the set.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    for (auto const& contig_id : contig_set) {

      if (not genome_ptr->getMap().contains(contig_id)) {

        auto contig_opt = genome_ptr->getCreateContig(contig_id);
        if (not contig_opt) {

          ExecEnv::log().error("PopulationDB::squareContigs; Unable to add contig: {} to genome: {}", contig_id, genome_id);

        }

      }

    }

  }

  return contig_set.size();

}


bool kgl::PopulationDB::addVariant( const std::shared_ptr<const Variant>& variant_ptr,
                                  const std::vector<GenomeId_t>& genome_vector) {

  bool result = true;
  for (auto& genome : genome_vector) {

    auto genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("PopulationDB::addVariant; Could not add/create genome: {}", genome);
      result = false;
      continue;

    }

    if (not genome_opt.value()->addVariant(variant_ptr)) {

      ExecEnv::log().error("PopulationDB::addVariant; Could not add variant to genome: {}", genome);
      result = false;
      continue;

    }

  }

  return result;

}


// Unconditional Merge.
size_t kgl::PopulationDB::mergePopulation(const std::shared_ptr<const PopulationDB>& merge_population) {

  size_t variant_count = 0;
  for (const auto& [genome, genome_ptr] : merge_population->getMap()) {

    auto genome_opt = getCreateGenome(genome);
    if (not genome_opt) {

      ExecEnv::log().error("PopulationDB::addUniqueVariant; Could not add/create genome: {}", genome);
      continue;
    }

    variant_count += genome_opt.value()->mergeGenome(genome_ptr);

  }

  return variant_count;

}


// Ensures that all variants are correctly specified.
std::pair<size_t, size_t> kgl::PopulationDB::validate(const std::shared_ptr<const GenomeReference>& genome_db) const {

  WorkflowThreads thread_pool(WorkflowThreads::defaultThreads());
  std::vector<std::future<std::pair<size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : getMap()) {

    // function, object_ptr, arg1
    std::future<std::pair<size_t, size_t>> future = thread_pool.enqueueFuture(&GenomeDB::validate,
                                                                              genome_ptr,
                                                                              genome_db);
    future_vector.push_back(std::move(future));

  }

  // Check the results of the validation.
  std::pair<size_t, size_t> population_count{0, 0};
  for (auto& future : future_vector) {

    std::pair<size_t, size_t> genome_count = future.get();

    if (genome_count.first != genome_count.second) {

      ExecEnv::log().warn("PopulationDB::validate(), Population: {} Failed to Validate, Total filter: {}, Validated: {}",
                          populationId(), genome_count.first, genome_count.second);

    }

    population_count.first += genome_count.first;
    population_count.second += genome_count.second;

  }

  return population_count;

}


std::shared_ptr<kgl::GenomeDB> kgl::PopulationDB::compressPopulation() const {

  std::shared_ptr<GenomeDB> compressedGenome(std::make_shared<GenomeDB>("Compressed"));

  if (not processAll(*compressedGenome, &GenomeDB::addVariant)) {

    ExecEnv::log().error("PopulationDB::compressPopulation(); problem compressing population: {}", populationId());
  }

  return  compressedGenome;

}


std::shared_ptr<kgl::PopulationDB> kgl::PopulationDB::uniqueUnphasedGenome() const {

  // Local class to perform the multi-threaded construction of a population with 1 genome containing
  // all the unique variants.
  class UniqueVariant {

  public:

    UniqueVariant(const PopulationDB& multi_genome) {

      compressed_population_ptr_ = std::make_shared<PopulationDB>(multi_genome.populationId() + "_Compressed", multi_genome.dataSource());

      auto unphased_genome_opt = compressed_population_ptr_->getCreateGenome("UniqueCompressed");

      if (not unphased_genome_opt) {

        ExecEnv::log().critical("PopulationDB::UniqueUnphased(); problem creating unique unphased genome with population: {}", multi_genome.populationId());

      }

      unphased_genome_ = unphased_genome_opt.value();

    }

    bool addUniqueUnphasedVariant(std::shared_ptr<const Variant> variant_ptr) {

      auto unphased_hash = variant_ptr->HGVS(); // Create a unique HGSV hash, phasing excluded.
      {
        // Acquire the mutex.
        std::scoped_lock lock(map_mutex_);

        // Check if the variant is already in the map.
        auto find_result = variant_map_.find(unphased_hash);
        if (find_result != variant_map_.end()) {

          // If already present, just return.
          return true;

        } else {

          // If not present, then add to the map.
          auto [insert_iter, result] = variant_map_.try_emplace(std::move(unphased_hash), variant_ptr);

          if (not result) {

            ExecEnv::log().error("PopulationDB::uniqueUnphasedGenome, cannot add duplicate variant hash: {}", variant_ptr->HGVS());
            return false;

          }

        }

      }
      // The variant was not found in the map so add to the unique population.
      bool add_result = unphased_genome_->addVariant(variant_ptr);
      if (not add_result) {

        ExecEnv::log().error("PopulationDB::uniqueUnphasedGenome, cannot add variant hash: {} to population", variant_ptr->HGVS());
        return false;

      }

      return true;

    }

    std::shared_ptr<PopulationDB> compressed_population_ptr_;
    std::shared_ptr<GenomeDB> unphased_genome_;
    // Implemented as a hash map for a bit of extra speed.
    std::unordered_map<std::string, std::shared_ptr<const Variant>> variant_map_;
    // Mutex to lock the map structure for safe multiple thread access.
    std::mutex map_mutex_;

  }; // End of local object definition.

  // Create an instance of the local object.
  std::shared_ptr<UniqueVariant> unique_variant_ptr(std::make_shared<UniqueVariant>(*this));

  // Using the local object and multi-threading, process all variants held in the population
  if (not processAll_MT(unique_variant_ptr, &UniqueVariant::addUniqueUnphasedVariant)) {

    ExecEnv::log().error("PopulationDB::UniqueUnphased(); problem creating unique unphased genome with population: {}", populationId());

  }

  // ReturnType the unique variant population.
  return  unique_variant_ptr->compressed_population_ptr_;

}


std::optional<std::shared_ptr<const kgl::InfoEvidenceHeader>> kgl::PopulationDB::getVCFInfoEvidenceHeader() const {

  for (auto const& genome : getMap()) {

    for (auto const& contig : genome.second->getMap()) {

      for (auto const& variant_vector : contig.second->getMap()) {

        for (auto const& variant_ptr : variant_vector.second->getVariantArray()) {

          if (variant_ptr->evidence().infoData()) {

            auto const& info_data = *(variant_ptr->evidence().infoData().value());
            return info_data.evidenceHeader();

          }

        }

      }

    }

  }

  return std::nullopt;

}

