//
// Created by kellerberrin on 6/07/23.
//

#include "kga_analysis_lib_seqmutation.h"
#include "kgl_variant_filter_db_offset.h"
#include "kgl_mutation_variant_filter.h"
#include "kel_workflow_threads.h"
#include "kgl_mutation_sequence.h"
#include "kgl_mutation_transcript.h"

#include <ranges>

namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kga::MutateGenes::mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr) {


  // Get the active contigs in this population.
  auto contig_map = population_ptr->contigCount();

  for (auto const& [contig_id, variant_count] : contig_map) {

    if (variant_count > 0) {

      auto gene_vector = contigGenes(contig_id);
      ExecEnv::log().info("Mutating contig_ref_ptr: {}, gene count: {}", contig_id, gene_vector.size());

      for (auto const& gene_ptr : gene_vector) {

        auto transcription_array = GeneFeature::getTranscriptionSequences(gene_ptr);
        for (auto const& [transcript_id,  transcript_ptr] : transcription_array->getMap()) {

          auto transcript_record = mutateTranscript( gene_ptr, transcript_id, transcript_ptr, population_ptr, genome_ptr_);
          mutate_analysis_.addTranscriptRecord(transcript_record);

        } // For transcript

      } // For genes.

    } // variants > 0

  } // For contigs.

}


std::vector<std::shared_ptr<const kgl::GeneFeature>> kga::MutateGenes::contigGenes(const ContigId_t& contig_id) const {

  std::vector<std::shared_ptr<const GeneFeature>> gene_vector;

  auto contig_opt = genome_ptr_->getContigSequence(contig_id);
  if (not contig_opt) {

    ExecEnv::log().warn("Requested contig: {} not found in reference genome: {}", contig_id, genome_ptr_->genomeId());
    return gene_vector;

  }
  auto const& contig_ref_ptr = contig_opt.value();
  for (auto const& [gene_offset, gene_ptr] : contig_ref_ptr->getGeneMap()) {

      gene_vector.push_back(gene_ptr);

  }

  return gene_vector;

}

// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
kga::TranscriptMutateRecord
kga::MutateGenes::mutateTranscript( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                    const FeatureIdent_t& transcript_id,
                                    const std::shared_ptr<const TranscriptionSequence>& transcript_ptr,
                                    const std::shared_ptr<const PopulationDB>& population_ptr,
                                    const std::shared_ptr<const GenomeReference>& reference_genome_ptr) {



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
  return transcript_record;

}


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
kga::MutateStats kga::MutateGenes::mutateGenomes( const std::shared_ptr<const GeneFeature>& gene_ptr,
                                                  const FeatureIdent_t& transcript_id,
                                                  const std::shared_ptr<const PopulationDB>& gene_population_ptr,
                                                  const std::shared_ptr<const GenomeReference>& reference_genome_ptr) {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  // The tuple is variant map, total modifying variants, multiple variants (more than 1 modifying variant per offset).
  using FutureType = std::future<std::pair<SequenceStats, bool>>;
  std::vector<std::pair<FutureType, std::string>> future_vector;


  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<std::pair<SequenceStats, bool>> future = thread_pool.enqueueFuture(&MutateGenes::genomeTranscriptMutation,
                                                                                this,
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
                                      stats.filter_statistics_.total_interval_variants_,
                                      stats.filter_statistics_.upstream_deleted_,
                                      stats.modified_sequence_,
                                      stats.modified_amino_size_,
                                      stats.original_sequence_,
                                      stats.original_amino_size_);

    // Calculate summary statistics.
    mutate_stats.total_variants_ += stats.filter_statistics_.total_interval_variants_;
    mutate_stats.total_snp_variants_ += stats.filter_statistics_.total_snp_variants_;
    mutate_stats.total_frameshift_variants_ += stats.filter_statistics_.total_frame_shift_;

    if (stats.filter_statistics_.non_unique_count_ > 0) {

      mutate_stats.duplicate_variants_ += stats.filter_statistics_.non_unique_count_;
      ++mutate_stats.duplicate_genomes_;

    }
    if (stats.filter_statistics_.total_interval_variants_ > 0) {

      ++mutate_stats.mutant_genomes_;

    }
    if (stats.filter_statistics_.upstream_deleted_ > 0) {

      mutate_stats.upstream_delete_variants_ += stats.filter_statistics_.upstream_deleted_;
      ++mutate_stats.upstream_delete_genomes_;

    }

    mutate_stats.modified_validity_.updateValidity(stats.modified_sequence_, stats.modified_amino_size_);
    mutate_stats.original_validity_.updateValidity(stats.original_sequence_, stats.original_amino_size_);

  } // For all futures (genomes).

  return mutate_stats;

}


std::pair<kga::SequenceStats, bool> kga::MutateGenes::genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_db_ptr,
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
  const auto& contig_ref_ptr = contig_ref_opt.value();

  // Get the gene gene_interval
  auto transcript_opt = contig_ref_ptr->getTranscription(gene_ptr->id(), transcript_id);
  if (not transcript_opt) {

    ExecEnv::log().warn("Transcript: {} not found, Gene: {}, Reference Genome: {}",
                        transcript_id, gene_ptr->id(), gene_contig_id);

    return {{}, false};

  }
  const auto& transcript_ptr = transcript_opt.value();

//  const auto filtertype = SeqVariantFilterType::HIGHEST_FREQ_VARIANT;
  const SequenceTranscript modified_transcript(contig_db_ptr, transcript_ptr, filtertype_);
  if (not modified_transcript.sequenceStatus()) {

    ExecEnv::log().warn("Unable to generate modified sequence for Transcript: {}, Gene: {}, Reference Genome: {}",
                        transcript_id, gene_ptr->id(), gene_contig_id);

    return {{}, false};

  }

  auto modified_opt = modified_transcript.getModifiedAdjustedValidity();
  auto original_opt = modified_transcript.getOriginalValidity();

  if (not modified_opt or not original_opt) {

    ExecEnv::log().warn("Problem retrieving modified sequence for Transcript: {}, Gene: {}, Reference Genome: {}",
                        transcript_id, gene_ptr->id(), gene_contig_id);

    return {{}, false};

  }

  auto& [modified_coding, modified_validity, modified_amino_size] = modified_opt.value();
  auto& [original_coding, original_validity, original_amino_size] = original_opt.value();

  SequenceStats stats;
  stats.filter_statistics_ = modified_transcript.filterStatistics();
  stats.modified_sequence_ = modified_validity;
  stats.modified_amino_size_ = modified_amino_size;
  stats.original_sequence_ = original_validity;
  stats.original_amino_size_ = original_amino_size;

  return { stats, true};

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kga::MutateGenesReport::mutatePopulation(const std::shared_ptr<const PopulationDB>& population_ptr) {

  // Mutate all the relevant genes in the relevant contigs.
  ExecEnv::log().info("MutateGenesReport; Begin gene mutation");
  mutate_genes_.mutatePopulation(population_ptr);
  ExecEnv::log().info("MutateGenesReport; End gene mutation");

}

void kga::MutateGenesReport::printMutateReports() {

  // Output the mutation statistics.
  std::string mutation_file_name = std::string("MutationTranscript") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, report_directory_);
  mutate_genes_.mutateAnalysis().printMutationTranscript(mutation_file_name);

  std::string validity_file_name = std::string("MutationValidity") + std::string(VARIANT_COUNT_EXT_);
  validity_file_name = Utility::filePath(validity_file_name, report_directory_);
  mutate_genes_.mutateAnalysis().printMutationValidity(validity_file_name);

  mutation_file_name = std::string("MutationGenome") + std::string(VARIANT_COUNT_EXT_);
  mutation_file_name = Utility::filePath(mutation_file_name, report_directory_);
  mutate_genes_.mutateAnalysis().printGenomeContig(mutation_file_name);

}
