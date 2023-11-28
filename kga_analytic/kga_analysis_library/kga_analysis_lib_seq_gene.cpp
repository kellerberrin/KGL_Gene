//
// Created by kellerberrin on 17/11/23.
//

#include "kel_workflow_threads.h"
#include "kga_analysis_lib_seq_gene.h"
#include "kgl_mutation_transcript.h"

#include <fstream>

namespace kel = kellerberrin;
namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
void kga::AnalysisTranscriptSequence::performTranscriptAnalysis( const std::shared_ptr<const PopulationDB>& gene_population_ptr) {

  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  // The tuple is variant map, total modifying variants, multiple variants (more than 1 modifying variant per offset).
  using FutureType = std::future<GenomeRecordOpt>;
  std::vector<std::pair<FutureType, std::string>> future_vector;


  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<GenomeRecordOpt> future = thread_pool.enqueueFuture(&AnalysisTranscriptSequence::genomeTranscriptMutation, this, genome_ptr);
    future_vector.push_back({std::move(future), genome_id});

  }

  // Wait for all the threads to return.
  for (auto & [future, genome_id] : future_vector) {

    auto record_opt  = future.get();
    if (not record_opt) {

      ExecEnv::log().warn("AnalysisTranscriptSequence; problem analyzing genome: {}, transcript: {}",
                          genome_id, transcript_ptr_->getParent()->id());
      // Skip the statistics
      continue;

    }
    auto const& record = record_opt.value();
    auto& [modified_coding, modified_validity, modified_amino_size] = record;
    auto string_view = modified_coding.getStringView();
    auto sequence_text = std::string(string_view);

    if (transcript_map_.contains(sequence_text)) {

      auto& [text, genomes] = *transcript_map_.find(sequence_text);
      genomes.genomes_.insert(genome_id);

    } else {

      TranscriptSequenceRecord transcript_record;
      transcript_record.genomes_.insert(genome_id);
      transcript_map_.try_emplace(sequence_text, transcript_record);

    }

  } // For all futures (genomes).



}


kga::GenomeRecordOpt kga::AnalysisTranscriptSequence::genomeTranscriptMutation(const std::shared_ptr<const GenomeDB>& genome_db_ptr) {

  // Reduce code noise.
  auto const& contig_ref_ptr = transcript_ptr_->getGene()->contig_ref_ptr();
  auto const& gene_contig_id = contig_ref_ptr->contigId();
  auto const& transcript_id = transcript_ptr_->getParent()->id();
  auto const& genome_id = genome_db_ptr->genomeId();

  // Use this to obtain the variant contig_ref_ptr.
  auto contig_db_opt = genome_db_ptr->getContig(gene_contig_id);
  if (not contig_db_opt) {

    ExecEnv::log().warn("Genome: {} does not contain variant contig: {} for transcript: {}", genome_id, gene_contig_id, transcript_id);
    return std::nullopt;

  }
  const auto& contig_db_ptr = contig_db_opt.value();

  const SequenceTranscript modified_transcript(contig_db_ptr, transcript_ptr_, FILTERTYPE_);
  if (not modified_transcript.sequenceStatus()) {

    ExecEnv::log().warn("Unable to generate modified sequence for Transcript: {}, Genome: {}", transcript_id, genome_id);
    return std::nullopt;

  }

  auto modified_opt = modified_transcript.getModifiedAdjustedValidity();
  auto original_opt = modified_transcript.getOriginalValidity();

  if (not modified_opt or not original_opt) {

    ExecEnv::log().warn("Problem retrieving modified sequence for Transcript: {}, Genome: {}", transcript_id, genome_id);
    return std::nullopt;

  }

//  auto& [modified_coding, modified_validity, modified_amino_size] = modified_opt.value();
//  auto& [original_coding, original_validity, original_amino_size] = original_opt.value();

  return modified_opt;

}


void kga::AnalysisTranscriptSequence::printReport(const std::string& report_directory) const {

  auto const& transcript_id = transcript_ptr_->getParent()->id();
  auto file_path = REPORT_PREFIX_ + "_" + transcript_id + REPORT_EXT_;
  file_path = kel::Utility::filePath(file_path, report_directory);

  std::ofstream report_stream(file_path);
  if (not report_stream.good()) {

    ExecEnv::log().warn("AnalysisTranscriptSequence; Could not open file: {}", file_path);
    return;

  }

  report_stream << "Genomes" << REPORT_FIELD_
                << "Text" << '\n';


  for (auto const& [text, record] : transcript_map_) {

    report_stream << record.genomes_.size() << REPORT_FIELD_
                  << text << '\n';

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kga::AnalysisTranscriptFamily::createAnalysisVector(const GeneVector& gene_vector) {

  // For all genes in the vector
  for (auto const& gene_ptr : gene_vector) {

    auto transcription_array_ptr = GeneFeature::getTranscriptionSequences(gene_ptr);
    auto const& transcription_map = transcription_array_ptr->getMap();

    // For all the transcripts of the gene.
    for (auto const& [transcript_id, transcript_ptr] : transcription_map) {

      analysis_vector_.emplace_back(transcript_ptr);

    }

  }

}

void kga::AnalysisTranscriptFamily::performFamilyAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr) {

  for (auto& transcript_analysis : analysis_vector_) {

    transcript_analysis.performTranscriptAnalysis(population_ptr);

  }

}


void kga::AnalysisTranscriptFamily::printAllReports(const std::string& analysis_directory,
                                                    const std::string& analysis_sub_directory) const {

  auto report_directory = Utility::appendPath(analysis_sub_directory, analysis_directory);

  // Recreate the report subdirectory.
  if (not Utility::directoryRenew(report_directory)) {

    ExecEnv::log().warn("AnalysisTranscriptFamily; cannot re-create report directory: {}", report_directory);

  }

  for (auto& transcript_analysis : analysis_vector_) {

    transcript_analysis.printReport(report_directory);

  }

}
