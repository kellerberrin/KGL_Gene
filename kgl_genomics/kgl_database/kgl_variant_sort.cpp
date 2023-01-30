//
// Created by kellerberrin on 27/6/21.
//

#include "kgl_variant_factory_vcf_evidence_analysis.h"
#include "kgl_variant_sort.h"
#include "../../kel_thread/kel_workflow_threads.h"


namespace kgl = kellerberrin::genome;




std::shared_ptr<kgl::EnsemblIndexMap> kgl::VariantSort::ensemblIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr) {

  std::shared_ptr<EnsemblIndexMap> ensembl_indexed_variants_ptr(std::make_shared<EnsemblIndexMap>());

  std::vector<std::string> empty_gene_list;

  ensemblAddIndex(unphased_population_ptr, empty_gene_list, ensembl_indexed_variants_ptr);

  return ensembl_indexed_variants_ptr;

}


// Index variants by the Ensembl gene code in the vep field.
// Only add the variants of the genes specified in the gene list.
// An empty list adds all variants.
void kgl::VariantSort::ensemblAddIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                          const std::vector<std::string>& ensembl_gene_list,
                                          std::shared_ptr<EnsemblIndexMap>& index_map_ptr) {

  // Local object performs the indexing.
  class IndexMap {

  public:

    IndexMap( std::shared_ptr<EnsemblIndexMap>& index_map_ptr,
              const std::vector<std::string>& ensembl_gene_list)
        : ensembl_indexed_variants_ptr_(index_map_ptr) {

      for (auto const& ensembl_gene_id : ensembl_gene_list) {

        ensembl_gene_set_.insert(ensembl_gene_id);

      }

    }
    ~IndexMap() = default;

    bool ensemblIndex(const std::shared_ptr<const Variant>& variant) {

      if (not initialized_) {

        field_index_ = InfoEvidenceAnalysis::getVepIndexes(*variant, std::vector<std::string>{VEP_ENSEMBL_FIELD_});
        initialized_ = true;

      }

      auto field_vector = InfoEvidenceAnalysis::getVepData(*variant, field_index_);

      // Only unique gene idents.
      for (auto const& field : field_vector) {

        // Only 1 field in the map.
        if (not field.empty()) {

          const auto& [field_ident, field_value] = *field.begin();
          if (not field_value.empty()) {

            unique_ident_.insert(field_value);

          }

        }

      }

      for (auto const& ident : unique_ident_) {

        if (not ident.empty()) {

          if (ensembl_gene_set_.empty() or ensembl_gene_set_.contains(ident)) {

            ensembl_indexed_variants_ptr_->emplace(ident, variant);

          }

        }

      }

      unique_ident_.clear();

      return true;

    }

  private:

    std::set<std::string> ensembl_gene_set_;
    std::shared_ptr<EnsemblIndexMap> ensembl_indexed_variants_ptr_;
    std::set<std::string> unique_ident_;
    VepIndexVector field_index_;
    bool initialized_{false};


  }; // IndexMap object.

  IndexMap index_map(index_map_ptr, ensembl_gene_list);

  unphased_population_ptr->processAll(index_map, &IndexMap::ensemblIndex);

}


size_t kgl::VariantSort::nonEnsemblIdentifiers(const EnsemblIndexMap& index_map) {

  size_t non_ensembl_identifiers{0};

  for (auto const& [ident, variant_ptr] : index_map) {

    if (ident.find(ENSEMBL_PREFIX_) == std::string::npos) {

      ++non_ensembl_identifiers;

    }

  }

  return non_ensembl_identifiers;

}


// Index by variant Id.
std::shared_ptr<kgl::VariantIdIndexMap> kgl::VariantSort::variantIdIndex(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Local object performs the indexing.
  class IndexMap {

  public:

    IndexMap() : indexed_variants_ptr_(std::make_shared<VariantIdIndexMap>()) {}
    ~IndexMap() = default;

    bool variantIDIndex(const std::shared_ptr<const Variant>& variant_ptr) {

      if (not variant_ptr->identifier().empty()) {

        indexed_variants_ptr_->emplace(variant_ptr->identifier(), variant_ptr);

      }

      return true;

    }

    std::shared_ptr<VariantIdIndexMap> indexedVariants() { return indexed_variants_ptr_; }

  private:

    std::shared_ptr<VariantIdIndexMap> indexed_variants_ptr_;

  }; // IndexMap object.

  IndexMap index_map;

  population_ptr->processAll(index_map, &IndexMap::variantIDIndex);

  return index_map.indexedVariants();

}


// Index by Genome and then by variant Id.
std::shared_ptr<kgl::VariantGenomeIndexMap> kgl::VariantSort::variantGenomeIndex(const std::shared_ptr<const PopulationDB>& population_ptr) {


  // Local object performs the indexing.
  class IndexMap {

  public:

    IndexMap() : indexed_variants_ptr_(std::make_shared<VariantIdIndexMap>()) {}
    ~IndexMap() = default;

    bool variantIDIndex(const std::shared_ptr<const Variant>& variant_ptr) {

      if (not variant_ptr->identifier().empty()) {

        indexed_variants_ptr_->emplace(variant_ptr->identifier(), variant_ptr);

      }

      return true;

    }

    std::shared_ptr<VariantIdIndexMap> indexedVariants() { return indexed_variants_ptr_; }

  private:

    std::shared_ptr<VariantIdIndexMap> indexed_variants_ptr_;

  }; // IndexMap object.

  std::shared_ptr<VariantGenomeIndexMap> genome_index_map(std::make_shared<VariantGenomeIndexMap>());

  for (auto const& [genome_id, genome_ptr] : population_ptr->getMap()) {

    IndexMap index_map;

    genome_ptr->processAll(index_map, &IndexMap::variantIDIndex);

    auto [iter, result] = genome_index_map->try_emplace(genome_id, index_map.indexedVariants());
    if (not result) {

      ExecEnv::log().error("VariantSort::variantGenomeIndex; Unable to (duplicate) add genome: {}", genome_id);

    }

  }

  return genome_index_map;

}



// Index by Genome and then by variant Id, multi-threaded.
std::shared_ptr<kgl::VariantGenomeIndexMap> kgl::VariantSort::variantGenomeIndexMT(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Local object performs the indexing.
  class IndexMap {

  public:

    IndexMap() : indexed_variants_ptr_(std::make_shared<VariantIdIndexMap>()) {}
    ~IndexMap() = default;

    bool variantIDIndex(const std::shared_ptr<const Variant>& variant_ptr) {

      if (not variant_ptr->identifier().empty()) {

        indexed_variants_ptr_->emplace(variant_ptr->identifier(), variant_ptr);

      }

      return true;

    }

    std::shared_ptr<VariantIdIndexMap> indexedVariants() { return indexed_variants_ptr_; }

    static std::shared_ptr<VariantIdIndexMap> indexGenome(const std::shared_ptr<const GenomeDB> genome_ptr) {

      IndexMap index_map;
      genome_ptr->processAll(index_map, &IndexMap::variantIDIndex);
      return index_map.indexedVariants();

    }

  private:

    std::shared_ptr<VariantIdIndexMap> indexed_variants_ptr_;

  }; // IndexMap object.


  // Thread count strategy
  size_t thread_count = WorkflowThreads::defaultThreads(population_ptr->getMap().size());
  // Fire-up the threads.
  WorkflowThreads thread_pool(thread_count);

  using GenomeIndexFuture = std::future<std::shared_ptr<VariantIdIndexMap>>;
  std::vector<std::pair<std::string, GenomeIndexFuture>> future_vector;

  for (auto const&[genome_id, genome_ptr] : population_ptr->getMap()) {

    GenomeIndexFuture future = thread_pool.enqueueTask(&IndexMap::indexGenome, genome_ptr);
    future_vector.emplace_back(genome_id, std::move(future));

  }

  std::shared_ptr<VariantGenomeIndexMap> genome_index_map(std::make_shared<VariantGenomeIndexMap>());
  // Retrieve the thread results into the map above.
  for (auto &future : future_vector) {

    auto& [genome_id, genome_future] = future;
    auto indexed_variant_ptr = genome_future.get();

    auto [iter, result] = genome_index_map->try_emplace(genome_id, indexed_variant_ptr);
    if (not result) {

      ExecEnv::log().error("VariantSort::variantGenomeIndexMT; Unable to add (duplicate) genome: {}", genome_id);

    }

  }

  return genome_index_map;

}

