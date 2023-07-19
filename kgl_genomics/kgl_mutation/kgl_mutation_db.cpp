//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_mutation_db.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db.h"
#include "kgl_variant_filter_unique.h"

#include <ranges>

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Used to find multiple different variants at a particular Genome/Contig/Offset.
// Multiple variants indicate different minor alleles at the same offset.
// P.Falciparum is haploid at the blood stage, so this indicates complexity of infection
// or an issue with the generation of the VCF file.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::MutateGenes::mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr) {


  // Get the active contigs in this population.
  auto contig_map = population_ptr->contigCount();

  for (auto const& [contig_id, variant_count] : contig_map) {

    if (variant_count > 0) {

      auto gene_vector = contigGenes(contig_id);
      ExecEnv::log().info("MutateGenes::mutatePopulation; Mutating contig: {}, gene count: {}", contig_id, gene_vector.size());

      for (auto const& gene_ptr : gene_vector) {

        auto transcription_array = GeneFeature::getTranscriptionSequences(gene_ptr);
        for (auto const& [transcript_id,  transcript_ptr] : transcription_array->getMap()) {

          mutateTranscript( gene_ptr, transcript_id, transcript_ptr, population_ptr, genome_ptr_);

        } // For transcript

      } // For genes.

    } // variants > 0

  } // For contigs.

}


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

// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
void kgl::MutateGenes::mutateTranscript( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                         const FeatureIdent_t& transcript_id,
                                         const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                         const std::shared_ptr<const PopulationDB>& population_ptr,
                                         const std::shared_ptr<const GenomeReference>& reference_genome) const {

  // Remove homozygous variants
  std::shared_ptr<const PopulationDB> unique_population_ptr = population_ptr->viewFilter(UniqueUnphasedFilter());
  // Only variants for the gene and transcript.
  std::shared_ptr<const PopulationDB> gene_population_ptr = unique_population_ptr->viewFilter(FilterGeneTranscriptVariants(gene_ptr, transcript_id));
  // Convert to a population with canonical variants.
  std::shared_ptr<const PopulationDB> canonical_population_ptr = gene_population_ptr->canonicalPopulation();
  // Remove any multiple variants.
  std::shared_ptr<const PopulationDB> mut_population_ptr = canonical_population_ptr->viewFilter(FrequencyUniqueFilter());

  // Check the structure of the canonical gene population.
  auto const [variant_count, validated_variants] = mut_population_ptr->validate(reference_genome);
  if (variant_count != validated_variants) {

    ExecEnv::log().error("MutateGenes::mutateTranscript; Failed Reference Genome Validation, Canonical Population Variant: {}, Canonical Validated: {}",
                         variant_count, validated_variants);

  }

  size_t mut_active_genome{0};
  for (auto const& [genome_id, genome_ptr] : mut_population_ptr->getMap()) {

    if (genome_ptr->variantCount() > 0) {

      ++mut_active_genome;

    }

  }

  auto const [ total_variants,
              duplicate_variants,
              duplicate_genomes] = mutateGenomes(gene_ptr, transcript_id, mut_population_ptr);

  ExecEnv::log().info("MutateGenes::mutateTranscript; Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("MutateGenes::mutateTranscript; Total Variants: {}, Processed Variants: {}, Duplicate Variants: {}, Total Genomes: {}, Mutant Genomes: {} Duplicate Genomes: {}",
                      mut_population_ptr->variantCount(),
                      total_variants,
                      duplicate_variants,
                      population_ptr->getMap().size(),
                      mut_active_genome,
                      duplicate_genomes);

  // Multiple Offset/variant statistics foe each transcript.
  TranscriptMutateRecord transcript_record( gene_ptr,
                                            transcript_ptr,
                                            total_variants,
                                            duplicate_variants,
                                            mut_active_genome,
                                            duplicate_genomes,
                                            population_ptr->getMap().size());

  mutate_analysis_.addTranscriptRecord(transcript_record);

}


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
std::tuple<size_t, size_t, size_t> kgl::MutateGenes::mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                                    const FeatureIdent_t& transcript_id,
                                                                    const std::shared_ptr<const PopulationDB>& gene_population_ptr) const {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  // the tuple is bool, genome_id, genome_id, transcript_id,
  std::vector<std::future<std::tuple<std::string, size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::tuple<std::string, size_t, size_t>> future = thread_pool.enqueueFuture(&MutateGenes::threadMutation, genome_ptr, gene_ptr, transcript_id);
    future_vector.push_back(std::move(future));

  }


  size_t total_variants{0};
  size_t duplicate_variants{0};
  size_t duplicate_genomes{0};
  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    auto const [genome_id, genome_total, genome_duplicate] = future.get();

    mutate_analysis_.addGenomeRecords(GenomeContigMutate(genome_id, genome_total, genome_duplicate));

    total_variants += genome_total;
    if (genome_duplicate > 0) {

      duplicate_variants += genome_duplicate;
      ++duplicate_genomes;

    }

  }

  return {total_variants, duplicate_variants, duplicate_genomes};

}

// .first total variants, .second multiple (duplicate) variants per offset.
std::tuple<std::string, size_t, size_t> kgl::MutateGenes::threadMutation( std::shared_ptr<const GenomeDB> genome_ptr,
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
        // Duplicate variant at offset.

        ++multiple_variants_;

      }

      return true;

    }

    size_t multipleVariants() const { return multiple_variants_; }
    std::map<ContigOffset_t, std::shared_ptr<const Variant>> ordered_variant_map_;

  private:

    std::shared_ptr<const GenomeDB> genome_ptr_;
    std::shared_ptr<const GeneFeature> gene_ptr_;
    size_t multiple_variants_{0};

  };

  ProcessVariants process_object(genome_ptr, gene_ptr);
  genome_ptr->processAll(process_object, &ProcessVariants::processVariants);

  return {genome_ptr->genomeId(), genome_ptr->variantCount(), process_object.multipleVariants() };

}

