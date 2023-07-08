//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_analysis_PfEMP_mutation.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db.h"


namespace kgl = kellerberrin::genome;


void kgl::MutateGenes::initializeGeneContigMap(const std::shared_ptr<const GenomeReference>& genome_ptr) {

  for (auto const& [contig_id, contig_ptr] : genome_ptr->getMap()) {

    for (auto const& [gene_offset, gene_ptr] : contig_ptr->getGeneMap()) {

      auto [iter, result] = gene_contig_map_.try_emplace(gene_ptr, contig_id);
      if (not result) {

        ExecEnv::log().warn("MutateGenes::initializeGeneContigMap; Unable to insert duplicate gene: {}", gene_ptr->id());

      }

    }

  }

}


std::vector<std::shared_ptr<const kgl::GeneFeature>> kgl::MutateGenes::contigGenes(const ContigId_t& contig_id) const {

  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;

  for (auto const& [gene_ptr, gene_contig_id] : gene_contig_map_) {

    if (gene_contig_id == contig_id) {

      gene_vector.push_back(gene_ptr);

    }

  }

  return gene_vector;

}


void  kgl::MutateGenes::mutateTranscript( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                          const FeatureIdent_t& transcript_id,
                                          const std::shared_ptr<const PopulationDB>& population_ptr,
                                          const std::shared_ptr<const GenomeReference>& reference_genome) const {

  // Remove homozygous variants
  std::shared_ptr<const PopulationDB> unique_population_ptr = population_ptr->viewFilter(UniqueUnphasedFilter());
  // Only variants for the gene and transcript.
  std::shared_ptr<const PopulationDB> gene_population_ptr = unique_population_ptr->viewFilter(FilterGeneTranscriptVariants(gene_ptr, transcript_id));
  // Convert to a population with canonical variants.
  std::shared_ptr<const PopulationDB> mut_population_ptr = gene_population_ptr->canonicalPopulation();

  // Checked the structure of the canonical gene population.
  auto const [variant_count, validated_variants] = mut_population_ptr->validate(reference_genome);
  if (variant_count != validated_variants) {

    ExecEnv::log().error("MutateGenes::mutateTranscript; Failed Reference Genome Validation, Canonical Population Variant: {}, Canonical Validated: {}",
                         variant_count, validated_variants);

  }

  ExecEnv::log().info("MutateGenes::mutateTranscript; Gene Population Variants: {}, Canonical Population Variant: {}, Canonical Validated: {}",
                      gene_population_ptr->variantCount(), variant_count, validated_variants);

  size_t mut_active_genome{0};
  for (auto const& [genome_id, genome_ptr] : mut_population_ptr->getMap()) {

    if (genome_ptr->variantCount() > 0) {

      ++mut_active_genome;

    }

  }

  size_t gene_active_genome{0};
  for (auto const& [genome_id, genome_ptr] : mut_population_ptr->getMap()) {

    if (genome_ptr->variantCount() > 0) {

      ++gene_active_genome;

    }

  }



  mutateGenomes(gene_ptr, transcript_id, mut_population_ptr);

  ExecEnv::log().info("MutateGenes::mutateTranscript; Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("MutateGenes::mutateTranscript; Variants: {}, Total Genomes: {}, Canonical Mutant Genomes: {}, Mutant Genomes: {}",
                      gene_population_ptr->variantCount(), gene_population_ptr->getMap().size(), mut_active_genome, gene_active_genome);

}



// Multi-tasked filtering for large populations.
void kgl::MutateGenes::mutateGenomes(const std::shared_ptr<const GeneFeature>& gene_ptr,
                                     const FeatureIdent_t& transcript_id,
                                     const std::shared_ptr<const PopulationDB>& gene_population_ptr) const {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  std::vector<std::future<bool>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<bool> future = thread_pool.enqueueFuture(&MutateGenes::threadMutation, genome_ptr, std::ref(gene_ptr), std::ref(transcript_id));
    future_vector.push_back(std::move(future));

  }


  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    bool  result = future.get();

  }

}

bool kgl::MutateGenes::threadMutation( std::shared_ptr<const GenomeDB> genome_ptr,
                                       const std::shared_ptr<const GeneFeature>& gene_ptr,
                                       const FeatureIdent_t& transcript_id) {

  class ProcessVariants {

  public:

    ProcessVariants(std::shared_ptr<const GenomeDB> genome_ptr,
                    std::shared_ptr<const GeneFeature> gene_ptr) : genome_ptr_(std::move(genome_ptr)), gene_ptr_(std::move(gene_ptr)) {}
    ~ProcessVariants() = default;

    bool processVariants(const std::shared_ptr<const Variant>& variant_ptr) {

      auto [ref_seq, alt_seq, canonical_offset] = variant_ptr->canonicalSequences();

      auto [insert_iter, result] = ordered_variant_map_.try_emplace(canonical_offset, variant_ptr);
      if (not result) {

        ExecEnv::log().warn("MutateGenes::threadMutation; Unable to add Variant: {}, Duplicate offset: {}, Genome: {}, Gene: {}"
                            ,  canonical_offset, variant_ptr->HGVS(), genome_ptr_->genomeId(), gene_ptr_->id());

        if (ordered_variant_map_.contains(canonical_offset)) {

          auto find_iter = ordered_variant_map_.find(canonical_offset);
          auto const& [map_offset, map_variant_ptr] = *find_iter;

          ExecEnv::log().warn("MutateGenes::threadMutation; Duplicate map offset: {}, Duplicate variant: {}, Genome:{}, Gene: {}",
                              map_offset, map_variant_ptr->HGVS(), genome_ptr_->genomeId(), gene_ptr_->id());

        } else {

          ExecEnv::log().error("MutateGenes::threadMutation; Unexpected map error, ordered map size: {}", ordered_variant_map_.size());

        }

      }

      return true;

    }
  private:

    std::map<ContigOffset_t, std::shared_ptr<const Variant>> ordered_variant_map_;
    std::shared_ptr<const GenomeDB> genome_ptr_;
    std::shared_ptr<const GeneFeature> gene_ptr_;

  };

  ProcessVariants process_object(genome_ptr, gene_ptr);
  genome_ptr->processAll(process_object, &ProcessVariants::processVariants);

  return true;

}