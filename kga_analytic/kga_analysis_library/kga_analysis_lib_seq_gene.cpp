//
// Created by kellerberrin on 17/11/23.
//

#include "kel_workflow_threads.h"
#include "kga_analysis_lib_seq_gene.h"
#include "kgl_mutation_transcript.h"
#include "kgl_sequence_distance_impl.h"

#include <fstream>
#include <ranges>

namespace kel = kellerberrin;
namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


// Multi-tasked filtering for large populations.
// .first total variants across all genomes, .second multiple (duplicate) variants per offset for all genomes.
void kga::AnalysisTranscriptSequence::performTranscriptAnalysis( const std::shared_ptr<const PopulationDB>& gene_population_ptr) {

  auto const& transcipt_id = transcript_ptr_->getParent()->id();
  // All other filters are multi-threaded for each genome.
  // Calc how many threads required.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  using FutureType = std::future<GenomeRecordOpt>;
  std::vector<std::pair<FutureType, std::string>> future_vector;


  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<GenomeRecordOpt> future = thread_pool.enqueueFuture(&AnalysisTranscriptSequence::genomeTranscriptMutation, this, genome_ptr);
    future_vector.emplace_back(std::move(future), genome_id);

  }

  // Wait for all the threads to return.
  size_t genome_count{0};
  for (auto & [future, genome_id] : future_vector) {

    auto record_opt  = future.get();
    if (not record_opt) {

      ExecEnv::log().warn("AnalysisTranscriptSequence; problem analyzing genome: {}, transcript: {}",
                          genome_id, transcipt_id);
      // Skip the statistics
      continue;

    }
    auto const& record = record_opt.value();
    auto const& [modified_coding, modified_validity, modified_amino_size] = record.modified();
    auto const& [original_coding, original_validity, original_amino_size] = record.reference();
    auto const& transcript_id = record.transcript()->getParent()->id();

    auto string_view = modified_coding.getStringView();
    auto sequence_text = std::string(string_view);

    if (transcript_map_.contains(sequence_text)) {

      auto& [text, sequence_record] = *transcript_map_.find(sequence_text);
      sequence_record.genomes_.insert(genome_id);
      sequence_record.transcripts_.insert(transcript_id);
      ++genome_count;

    } else {

      TranscriptSequenceRecord transcript_record;
      double distance = LevenshteinGlobalCoding(modified_coding, original_coding);
      transcript_record.distance_ = distance;
      transcript_record.genomes_.insert(genome_id);
      transcript_record.transcripts_.insert(transcript_id);
      auto [insert_iter, result] = transcript_map_.try_emplace(sequence_text, transcript_record);
      if (not result) {

        ExecEnv::log().error("Unable to insert sequence record for genome: {}", genome_id);

      } else {

        ++genome_count;

      }

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

  TranscriptModifyRecord modify_record(genome_id, transcript_ptr_, std::move(original_opt.value()), std::move(modified_opt.value()));
  return modify_record;

}


void kga::AnalysisTranscriptSequence::printReport(const std::string& report_directory) const {

  auto const& transcript_id = transcript_ptr_->getParent()->id();
  printReport(transcript_id, report_directory);

}

void kga::AnalysisTranscriptSequence::printReport(const std::string& report_label, const std::string& report_directory) const {

  auto file_path = REPORT_PREFIX_ + "_" + report_label + REPORT_EXT_;
  file_path = kel::Utility::filePath(file_path, report_directory);

  std::ofstream report_stream(file_path);
  if (not report_stream.good()) {

    ExecEnv::log().warn("AnalysisTranscriptSequence; Could not open file: {}", file_path);
    return;

  }

  report_stream << "Genomes" << REPORT_FIELD_
                << "TranscriptCount" << REPORT_FIELD_
                << "Transcripts" << REPORT_FIELD_
                << "Distance" << REPORT_FIELD_
                << "Text" << '\n';


  for (auto const& [text, record] : transcript_map_) {

    std::string transcript_list;
    for (auto const& transcript_id : record.transcripts_) {

      if (not transcript_list.empty()) {

        transcript_list += REPORT_SUBFIELD_;

      }

      transcript_list += transcript_id;

    }

    report_stream << record.genomes_.size() << REPORT_FIELD_
                  << record.transcripts_.size() << REPORT_FIELD_
                  << transcript_list << REPORT_FIELD_
                  << record.distance_ << REPORT_FIELD_
                  << text << '\n';

  }

}

void kga::AnalysisTranscriptSequence::mergeMap(const TranscriptMap& transcript_map) {

  for (auto const& [sequence_text, merge_record] : transcript_map) {

    if (transcript_map_.contains(sequence_text)) {

      auto& [text, map_record] = *transcript_map_.find(sequence_text);
      map_record.genomes_.insert(merge_record.genomes_.begin(), merge_record.genomes_.end());
      map_record.transcripts_.insert(merge_record.transcripts_.begin(), merge_record.transcripts_.end());

    } else {

      auto [insert_iter, result] = transcript_map_.try_emplace(sequence_text, merge_record);
      if (not result) {

        ExecEnv::log().error("AnalysisTranscriptSequence::mergeMap; Unable to merge sequence record");

      }

    }

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

    // Only process the transcript if the population contains variants for the transcript contig.
    const auto& transcript = transcript_analysis.transcript();
    const auto& transcript_contig = transcript.contig()->contigId();
    auto contig_count_opt = population_ptr->contigCount(transcript_contig);
    if (not contig_count_opt) {

      ExecEnv::log().warn("Unexpected, population: {} does not contain the contig: {} for transcript: {}",
                          population_ptr->populationId(), transcript_contig, transcript.getParent()->id());
      continue;

    }
    const auto& contig_count = contig_count_opt.value();
    if (contig_count > 0) {

      transcript_analysis.performTranscriptAnalysis(population_ptr);

    }

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

  if (not analysis_vector_.empty()) {

    auto merged_analysis = generateTotal();
    merged_analysis.printReport("MergedTotal", report_directory);

  }

}

kga::AnalysisTranscriptSequence kga::AnalysisTranscriptFamily::generateTotal() const {

  // Copy the first view
  AnalysisTranscriptSequence merged_sequence(analysis_vector_.front());
  // Vector will preserve sort order, drop first view
  for (auto const& trans_sequence : std::ranges::drop_view{ analysis_vector_, 1}) {

    merged_sequence.mergeMap(trans_sequence.transcriptMap());

  }

  return merged_sequence;

}
