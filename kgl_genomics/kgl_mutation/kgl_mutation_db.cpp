//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_mutation_db.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_variant_filter_unique.h"
#include "kgl_mutation_offset.h"
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
  // The tuple is variant map, total modifying variants, multiple variants (more than 1 modifying variant per offset).
  std::vector<std::future<std::tuple<std::shared_ptr<const RegionVariantMap>, size_t, size_t>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::tuple<std::shared_ptr<const RegionVariantMap>, size_t, size_t>> future = thread_pool.enqueueFuture(&MutateGenes::genomeTranscriptMutation, genome_ptr, gene_ptr, transcript_id);
    future_vector.push_back(std::move(future));

  }

  size_t total_variants{0}; // Total variants.
  size_t duplicate_variants{0}; // More than 1 variant per offset
  size_t duplicate_genomes{0}; // Genomes containing more than 1 variant per offset.
  size_t mutant_genomes{0}; // Genomes with mutations.
  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    auto const [unique_transcript_ptr, variant_count, duplicates] = future.get();

    mutate_analysis_.addGenomeRecords(GenomeContigMutate(unique_transcript_ptr->genomeId(), variant_count, duplicates));

    total_variants += variant_count;
    if (duplicates > 0) {

      duplicate_variants += duplicates;
      ++duplicate_genomes;

    }
    if (variant_count > 0) {

      ++mutant_genomes;

    }

  }

  return {total_variants, mutant_genomes, duplicate_variants, duplicate_genomes};

}


std::tuple<std::shared_ptr<const kgl::RegionVariantMap>, size_t, size_t>
kgl::MutateGenes::genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_ptr,
                                           const std::shared_ptr<const GeneFeature>& gene_ptr,
                                           const FeatureIdent_t& transcript_id) {

  ContigId_t gene_contig = gene_ptr->contig()->contigId();
  auto contig_opt = genome_ptr->getContig(gene_contig);
  if (not contig_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; Genome: {} does not contain contig: {} for gene: {}",
                        genome_ptr->genomeId(), gene_contig, gene_ptr->id());

    auto unique_transcript_ptr = std::make_shared<const RegionVariantMap>(genome_ptr->genomeId(),
                                                                          "",
                                                                          gene_ptr->sequence().begin(),
                                                                          gene_ptr->sequence().end(),
                                                                          OffsetVariantMap());
    return { unique_transcript_ptr, 0, 0 };

  }
  auto contig_ptr = contig_opt.value();

  GeneIntervalStructure gene_interval(gene_ptr);

  auto interval = gene_interval.geneInterval();

  auto [interval_map, non_unique_count] = MutationOffset::getCanonicalVariants(contig_ptr, interval.lower(), interval.upper());

  size_t map_size = interval_map.size() + non_unique_count;

  auto unique_transcript_ptr = std::make_shared<const RegionVariantMap>(genome_ptr->genomeId(),
                                                                        contig_ptr->contigId(),
                                                                        gene_ptr->sequence().begin(),
                                                                        gene_ptr->sequence().end(),
                                                                        std::move(interval_map));

  return {unique_transcript_ptr, map_size, non_unique_count };

}


