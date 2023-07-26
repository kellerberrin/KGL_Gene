//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_mutation_db.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db.h"
#include "kgl_variant_filter_unique.h"
#include "kel_workflow_threads.h"


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


  // Remove homozygous variants
  std::shared_ptr<const PopulationDB> unique_population_ptr = population_ptr->viewFilter(UniqueUnphasedFilter());

  // Get the active contigs in this population.
  auto contig_map = unique_population_ptr->contigCount();

  for (auto const& [contig_id, variant_count] : contig_map) {

    if (variant_count > 0) {

      auto gene_vector = contigGenes(contig_id);
      ExecEnv::log().info("MutateGenes::mutatePopulation; Mutating contig: {}, gene count: {}", contig_id, gene_vector.size());

      for (auto const& gene_ptr : gene_vector) {

        auto transcription_array = GeneFeature::getTranscriptionSequences(gene_ptr);
        for (auto const& [transcript_id,  transcript_ptr] : transcription_array->getMap()) {

          mutateTranscript( gene_ptr, transcript_id, transcript_ptr, unique_population_ptr, genome_ptr_);

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



  auto const [ total_variants,
               mut_active_genome,
              duplicate_variants,
              duplicate_genomes] = mutateGenomes(gene_ptr, transcript_id, population_ptr);

  ExecEnv::log().info("MutateGenes::mutateTranscript; Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("MutateGenes::mutateTranscript; Processed Variants: {}, Duplicate Variants: {}, Total Genomes: {}, Mutant Genomes: {} Duplicate Genomes: {}",
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
std::tuple<size_t, size_t, size_t, size_t> kgl::MutateGenes::mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
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

    std::future<std::tuple<std::string, size_t, size_t>> future = thread_pool.enqueueFuture(&MutateGenes::genomeTranscriptMutation, genome_ptr, gene_ptr, transcript_id);
    future_vector.push_back(std::move(future));

  }

  size_t total_variants{0};
  size_t duplicate_variants{0};
  size_t duplicate_genomes{0};
  size_t mutant_genomes{0};
  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    auto const [genome_id, genome_total, genome_duplicate] = future.get();

    mutate_analysis_.addGenomeRecords(GenomeContigMutate(genome_id, genome_total, genome_duplicate));

    total_variants += genome_total;
    if (genome_duplicate > 0) {

      duplicate_variants += genome_duplicate;
      ++duplicate_genomes;

    }
    if (genome_total > 0) {

      ++mutant_genomes;

    }

  }

  return {total_variants, mutant_genomes, duplicate_variants, duplicate_genomes};

}

// .first total variants, .second multiple (duplicate) variants per offset.
std::tuple<std::string, size_t, size_t>
    kgl::MutateGenes::genomeTranscriptMutation( const std::shared_ptr<const GenomeDB>& genome_ptr,
                                                const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                const FeatureIdent_t& transcript_id) {

  class ProcessVariants {

  public:

    ProcessVariants(std::shared_ptr<const GeneFeature> gene_ptr) : gene_ptr_(std::move(gene_ptr)) {}
    ~ProcessVariants() = default;

    bool processVariants(const std::shared_ptr<const Variant>& variant_ptr) {

      if (not variant_ptr->isCanonical()) {

        ExecEnv::log().error("MutateGenes::genomeTranscriptMutation; variant: {} is not canonical", variant_ptr->HGVS());

      }

      ContigOffset_t map_offset = variant_ptr->offset();
      // If an indel then insert into the next offset.
      if (not variant_ptr->isSNP()) {

        ++map_offset;

      }

      auto [insert_iter, result] = ordered_variant_map_.try_emplace(map_offset, variant_ptr);
      if (not result) {

        // Duplicate variant at offset.
        ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; duplicate variant: {}", variant_ptr->HGVS());
        ++multiple_variants_;

      }

      return true;

    }

    [[nodiscard]] size_t multipleVariants() const { return multiple_variants_; }
    OffsetVariantMap ordered_variant_map_;

  private:

    std::shared_ptr<const GeneFeature> gene_ptr_;
    size_t multiple_variants_{0};

  };

  ContigId_t gene_contig = gene_ptr->contig()->contigId();
  auto contig_opt = genome_ptr->getContig(gene_contig);
  if (not contig_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; Genome: {} does not contain contig: {} for gene: {}",
                        genome_ptr->genomeId(), gene_contig, gene_ptr->id());
    return { genome_ptr->genomeId(), 0, 0 };

  }


  // Only variants for the gene and transcript.
  std::shared_ptr<const ContigDB> gene_contig_ptr = contig_opt.value()->viewFilter(FilterGeneTranscriptVariants(gene_ptr, transcript_id));
  // Convert to a contig with canonical variants.
  std::shared_ptr<const ContigDB> canonical_contig_ptr = gene_contig_ptr->canonicalContig();
  // Remove any multiple mutating variants.
  std::shared_ptr<const ContigDB> unique_contig_ptr = canonical_contig_ptr->viewFilter(FrequencyUniqueFilter());

  size_t canonical_count = canonical_contig_ptr->variantCount();
  size_t duplicate_variants = canonical_count - unique_contig_ptr->variantCount();

  ProcessVariants process_object(gene_ptr);
  unique_contig_ptr->processAll(process_object, &ProcessVariants::processVariants);

  return {genome_ptr->genomeId(), canonical_count, duplicate_variants };

}

