//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_mutation_db.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_mutation_offset.h"
#include "kel_workflow_threads.h"
#include "kgl_mutation_sequence.h"
#include "kgl_mutation_aggregation.h"

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
                                         const std::shared_ptr<const GenomeReference>& reference_genome_ptr) const {



  auto const mutate_stats = mutateGenomes(gene_ptr, transcript_id, population_ptr, reference_genome_ptr);

  ExecEnv::log().info("MutateGenes::mutateTranscript; Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("MutateGenes::mutateTranscript; Processed Variants: {}, Duplicate Variants: {}, Total Genomes: {}, Mutant Genomes: {} Duplicate Genomes: {}, Upstream Deleted: {}",
                      mutate_stats.total_variants_,
                      mutate_stats.duplicate_variants_,
                      population_ptr->getMap().size(),
                      mutate_stats.mutant_genomes_,
                      mutate_stats.duplicate_genomes_,
                      mutate_stats.upstream_delete_variants_);

  // Multiple Offset/variant statistics for each transcript.
  TranscriptMutateRecord transcript_record( gene_ptr, transcript_ptr, mutate_stats);
  mutate_analysis_.addTranscriptRecord(transcript_record);

}


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
kgl::MutateStats kgl::MutateGenes::mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                  const FeatureIdent_t& transcript_id,
                                                  const std::shared_ptr<const PopulationDB>& gene_population_ptr,
                                                  const std::shared_ptr<const GenomeReference>& reference_genome_ptr) const {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  // The tuple is variant map, total modifying variants, multiple variants (more than 1 modifying variant per offset).
  std::vector<std::future<std::optional<RegionReturn>>> future_vector;

  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::optional<RegionReturn>> future = thread_pool.enqueueFuture(&MutateGenes::genomeTranscriptMutation,
                                                                                genome_ptr,
                                                                                gene_ptr,
                                                                                transcript_id,
                                                                                reference_genome_ptr);
    future_vector.push_back(std::move(future));

  }

  MutateStats mutate_stats;
  // Wait for all the threads to return.
  for (auto& future : future_vector) {

    const std::optional<RegionReturn> region_return_opt  = future.get();

    if (not region_return_opt) {

      // Skip the statistics
      continue;

    }

    ++mutate_stats.total_genomes_;

    auto const& region_return = region_return_opt.value();

    // Collect statistics.
    GenomeContigMutate genome_contig(region_return.regionVariant()->genomeId(),
                                     region_return.variantCount(),
                                     region_return.duplicates());
    mutate_analysis_.addGenomeRecords(genome_contig);

    // Calculate summary statistics.
    mutate_stats.total_variants_ += region_return.variantCount();
    if (region_return.duplicates() > 0) {

      mutate_stats.duplicate_variants_ += region_return.duplicates();
      ++mutate_stats.duplicate_genomes_;

    }
    if (region_return.variantCount() > 0) {

      ++mutate_stats.mutant_genomes_;

    }
    if (region_return.upstreamDeleted() > 0) {

      mutate_stats.upstream_delete_variants_ += region_return.upstreamDeleted();
      ++mutate_stats.upstream_delete_genomes_;

    }

  }

  return mutate_stats;

}


std::optional<kgl::RegionReturn> kgl::MutateGenes::genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_ptr,
                                                                            const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                                            const FeatureIdent_t& transcript_id,
                                                                            const std::shared_ptr<const GenomeReference>& reference_genome_ptr) {

  // Get the gene contig id.
  const ContigId_t gene_contig_id = gene_ptr->contig()->contigId();

  // Use this to obtain the variant contig.
  auto contig_opt = genome_ptr->getContig(gene_contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; Genome: {} does not contain contig: {} for gene: {}",
                        genome_ptr->genomeId(), gene_contig_id, gene_ptr->id());

    return std::nullopt;

  }
  auto contig_ptr = contig_opt.value();

  // And the reference genome contig.
  auto contig_ref_opt = reference_genome_ptr->getContigSequence(gene_contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; Reference Genome: {} does not contain contig: {} for gene: {}",
                        reference_genome_ptr->genomeId(), gene_contig_id, gene_ptr->id());

    return std::nullopt;

  }
  auto contig_ref_ptr = contig_ref_opt.value();

  // Get the gene interval
  GeneIntervalStructure gene_interval(gene_ptr);
  auto interval = gene_interval.geneInterval();

  // Return filtered variants adjusted for duplicate variants and upstream deletes.
  auto [interval_map, non_unique_count, upstream_deleted] = MutationOffset::getCanonicalVariants(contig_ptr, interval.lower(), interval.upper());

  size_t map_size = interval_map.size() + non_unique_count;

  auto region_variant_ptr = std::make_shared<const RegionVariantMap>(genome_ptr->genomeId(),
                                                                     contig_ptr->contigId(),
                                                                     gene_ptr->sequence().begin(),
                                                                     gene_ptr->sequence().end(),
                                                                     std::move(interval_map));

  auto adjusted_offset_ptr = std::make_shared<AdjustedSequenceInterval>(region_variant_ptr->variantRegion());
  if (not adjusted_offset_ptr->processVariantMap(region_variant_ptr->variantMap())) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation;  transcript: {}, problem updating interval: {}, genome: {}, contig: {}",
                        transcript_id,
                        region_variant_ptr->variantRegion().toString(),
                        region_variant_ptr->genomeId(),
                        region_variant_ptr->contigId());
    return std::nullopt;

  }

  auto adjusted_sequence_ptr = std::make_shared<AdjustedSequence>( contig_ref_ptr,
                                                                   region_variant_ptr->variantRegion(),
                                                                   adjusted_offset_ptr->indelModifyMap());

  if (not adjusted_sequence_ptr->updateSequence()) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; transcript: {}, problem updating sequence interval: {} genome: {}, contig: {}",
                        transcript_id,
                        region_variant_ptr->variantRegion().toString(),
                        region_variant_ptr->genomeId(),
                        region_variant_ptr->contigId());
    return std::nullopt;

  }

  auto modified_sequence_opt = SequenceAggregation::getModifiedGene(*adjusted_sequence_ptr, gene_interval, transcript_id);
  if (not modified_sequence_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; problem calculating modified gene sequence for transcript : {}, gene: {}",
                        transcript_id, gene_ptr->id());

  }

  auto original_sequence_opt = SequenceAggregation::getOriginalGene(*adjusted_sequence_ptr, gene_interval, transcript_id);
  if (not original_sequence_opt) {

    ExecEnv::log().warn("MutateGenes::genomeTranscriptMutation; problem calculating original gene sequence for transcript : {}, gene: {}",
                        transcript_id, gene_ptr->id());

  }

  if (modified_sequence_opt and original_sequence_opt) {

    auto& modify_sequence = modified_sequence_opt.value();
    auto& original_sequence = original_sequence_opt.value();

    if (modify_sequence.length() == 0) {

      ExecEnv::log().info("***** MutateGenes::genomeTranscriptMutation; deleted gene: {}, transcript: {}, original size: {}",
                          gene_ptr->id(), transcript_id, original_sequence.length());

    }

  }


  RegionReturn region_return(region_variant_ptr, map_size, non_unique_count, upstream_deleted);

  return region_return;

}
