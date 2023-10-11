//
// Created by kellerberrin on 6/07/23.
//

#include "kgl_mutation_db.h"
#include "kgl_variant_filter_features.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_genome_seq/kgl_seq_variant_filter.h"
#include "kel_workflow_threads.h"
#include "kgl_mutation_sequence.h"
#include "kgl_genome_seq/kgl_seq_transcript.h"

#include <ranges>

namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::MutateGenes::mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr) {


  // Get the active contigs in this population.
  auto contig_map = population_ptr->contigCount();

  for (auto const& [contig_id, variant_count] : contig_map) {

    if (variant_count > 0) {

      auto gene_vector = contigGenes(contig_id);
      ExecEnv::log().info("Mutating contig_ref_ptr: {}, gene count: {}", contig_id, gene_vector.size());

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

        ExecEnv::log().warn("Unable to insert duplicate gene: {}", gene_ptr->id());

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

  ExecEnv::log().info("Filtered Gene: {}, Transcript: {}, Description: {}",
                      gene_ptr->id(), transcript_id, gene_ptr->descriptionText());
  ExecEnv::log().info("Processed Variants: {}, Duplicate Variants: {}, Total Genomes: {}, Mutant Genomes: {} Duplicate Genomes: {}, Upstream Deleted: {}",
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
  using FutureType = std::future<std::pair<kgl::SequenceStats, bool>>;
  std::vector<std::pair<FutureType, std::string>> future_vector;


  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::pair<kgl::SequenceStats, bool>> future = thread_pool.enqueueFuture(&MutateGenes::genomeTranscriptMutation,
                                                                                genome_ptr,
                                                                                gene_ptr,
                                                                                transcript_id,
                                                                                reference_genome_ptr);
    future_vector.push_back({std::move(future), genome_id});

  }

  MutateStats mutate_stats;
  // Wait for all the threads to return.
  for (auto & [future, genome_id] : future_vector) {

    auto const [stats, valid_return]  = future.get();
    if (not valid_return) {

      // Skip the statistics
      continue;

    }

    ++mutate_stats.total_genomes_;

    // Collect statistics.
    mutate_analysis_.addGenomeRecords(genome_id,
                                      stats.map_size_,
                                      stats.upstream_deleted_,
                                      stats.modified_sequence_,
                                      stats.original_sequence_);

    // Calculate summary statistics.
    mutate_stats.total_variants_ += stats.map_size_;
    if (stats.non_unique_count_ > 0) {

      mutate_stats.duplicate_variants_ += stats.non_unique_count_;
      ++mutate_stats.duplicate_genomes_;

    }
    if (stats.map_size_ > 0) {

      ++mutate_stats.mutant_genomes_;

    }
    if (stats.upstream_deleted_ > 0) {

      mutate_stats.upstream_delete_variants_ += stats.upstream_deleted_;
      ++mutate_stats.upstream_delete_genomes_;

    }

    mutate_stats.modified_validity_.updateValidity(stats.modified_sequence_);
    mutate_stats.original_validity_.updateValidity(stats.original_sequence_);

  }

  return mutate_stats;

}


std::pair<kgl::SequenceStats, bool> kgl::MutateGenes::genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_db_ptr,
                                                                               const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                                               const FeatureIdent_t& transcript_id,
                                                                               const std::shared_ptr<const GenomeReference>& genome_ref_ptr) {

  // Get the gene contig_ref_ptr id.
  const ContigId_t gene_contig_id = gene_ptr->contig_ref_ptr()->contigId();

  // Use this to obtain the variant contig_ref_ptr.
  auto contig_db_opt = genome_db_ptr->getContig(gene_contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Genome: {} does not contain variant contig: {} for gene: {}",
                        genome_db_ptr->genomeId(), gene_contig_id, gene_ptr->id());

    return {{}, false};

  }
  const auto& contig_db_ptr = contig_db_opt.value();

  // And the reference genome contig_ref_ptr.
  auto contig_ref_opt = genome_ref_ptr->getContigSequence(gene_contig_id);
  if (not contig_ref_opt) {

    ExecEnv::log().warn("Reference Genome: {} does not contain reference contig: {} for gene: {}",
                        genome_ref_ptr->genomeId(), gene_contig_id, gene_ptr->id());

    return {{}, false};

  }
  const auto contig_ref_ptr = contig_ref_opt.value();

  // Get the gene gene_interval
  auto transcript_opt = contig_ref_ptr->getTranscription(gene_ptr->id(), transcript_id);
  if (not transcript_opt) {

    ExecEnv::log().warn("Transcript: {} not found, Gene: {}, Reference Genome: {}",
                        transcript_id, gene_ptr->id(), gene_contig_id);

    return {{}, false};

  }
  const auto transcript_ptr = transcript_opt.value();

  SequenceTranscript modified_sequence(contig_db_ptr, transcript_ptr);
  if (not modified_sequence.sequenceStatus()) {

    ExecEnv::log().warn("Unable to generate modified seuquence for Transcript: {}, Gene: {}, Reference Genome: {}",
                        transcript_id, gene_ptr->id(), gene_contig_id);

    return {{}, false};

  }

  auto stats = modified_sequence.sequenceStatistics();

  auto modified_sequence_opt = modified_sequence.getModifiedGene();
  auto original_sequence_opt = modified_sequence.getOriginalGene();

  if (modified_sequence_opt and original_sequence_opt) {

    if (transcript_ptr->codingType() == TranscriptionSequenceType::PROTEIN) {

      auto& modified = modified_sequence_opt.value();
      auto modified_coding = modified.codingSequence(transcript_ptr->strand());
      stats.modified_sequence_ = gene_ptr->contig_ref_ptr()->checkValidCodingSequence(modified_coding);

      auto& original = original_sequence_opt.value();
      auto original_coding = original.codingSequence(transcript_ptr->strand());
      stats.original_sequence_ = gene_ptr->contig_ref_ptr()->checkValidCodingSequence(original_coding);

    } else {

      stats.modified_sequence_ = CodingSequenceValidity::NCRNA;
      stats.original_sequence_ = CodingSequenceValidity::NCRNA;

    }

  } else {

    ExecEnv::log().warn("Problem mutating Gene: {}, Transcript: {}, Genome: {}",
                        gene_ptr->id(), transcript_id, genome_db_ptr->genomeId());

  }

  return { stats, true};

}
