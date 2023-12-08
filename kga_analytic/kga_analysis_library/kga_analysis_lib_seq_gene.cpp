//
// Created by kellerberrin on 17/11/23.
//

#include "kel_workflow_threads.h"
#include "kga_analysis_lib_seq_gene.h"
#include "kgl_mutation_transcript.h"
#include "kgl_sequence_distance_impl.h"

#include <fstream>
#include <ranges>
#include <functional>

namespace kel = kellerberrin;
namespace kgl = kellerberrin::genome;
namespace kga = kellerberrin::genome::analysis;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kga::TranscriptSequenceRecord::TranscriptSequenceRecord(const std::shared_ptr<const TranscriptModifyRecord>& genome_modify_ptr)
: modified_sequence_(genome_modify_ptr->modified().getView()) {

  if (not addSequenceRecord(genome_modify_ptr)) {

    ExecEnv::log().warn("Unable to add transcript sequence record for genome: {}", genome_modify_ptr->genome());

  }

}


bool kga::TranscriptSequenceRecord::addSequenceRecord(const std::shared_ptr<const TranscriptModifyRecord>& genome_modify_ptr) {

  // Check the sequence views match.
  if (genome_modify_ptr->modified().getView() != modified_sequence_) {

    ExecEnv::log().warn("Cannot add transcript sequence record for genome: {}", genome_modify_ptr->genome());
    auto add_id = genome_modify_ptr->transcript()->getParent()->id();
    auto add_view = genome_modify_ptr->modified().getView();
    ExecEnv::log().warn("Sequence: {} expected size: {}, add sequence size: {}",
                        add_id,
                        modified_sequence_.length(),
                        add_view.length());
    ExecEnv::log().warn("Sequence: {} hash: {}, add sequence hash: {}",
                        add_id,
                        TranscriptSequenceRecord::generateSequenceLabel(modified_sequence_),
                        TranscriptSequenceRecord::generateSequenceLabel(add_view));
    return false;

  }

  // Add to the genome map.
  bool result{true};
  genomes_.emplace(genome_modify_ptr->genome(), genome_modify_ptr);

  // Add to the transcript map.
  auto const& transcript_id = genome_modify_ptr->transcript()->getParent()->id();
  if (transcripts_.contains(transcript_id)) {

    auto& [find_id, transcript_vector] = *transcripts_.find(transcript_id);
    transcript_vector.push_back(genome_modify_ptr);

  } else {

    auto [transcript_iter, transcript_result] = transcripts_.try_emplace(transcript_id, TransSeqVector{genome_modify_ptr});
    if (not transcript_result) {

      ExecEnv::log().warn("Cannot add transcript sequence record to genome map, genome: {}", genome_modify_ptr->genome());
      result = false;

    }

  }

  return result;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kga::AnalysisTranscriptSequence::performTranscriptAnalysis( const std::shared_ptr<const PopulationDB>& gene_population_ptr) {

  auto const& transcipt_id = transcript_ptr_->getParent()->id();
  // Multithreading for each genome.
  // Calc how many threads are required as a function of available HW processors.
  size_t thread_count = WorkflowThreads::defaultThreads(gene_population_ptr->getMap().size());
  WorkflowThreads thread_pool(thread_count);
  // A vector for futures.
  using FutureType = std::future<GenomeRecordOpt>;
  std::vector<FutureType> future_vector;


  // Queue a thread for each genome.
  for (auto const& [genome_id, genome_ptr] : gene_population_ptr->getMap()) {

    std::future<GenomeRecordOpt> future = thread_pool.enqueueFuture(&AnalysisTranscriptSequence::genomeTranscriptMutation, this, genome_ptr);
    future_vector.emplace_back(std::move(future));

  }

  // Process the queued futures as they become available.
  for (auto & future : future_vector) {

    auto record_ptr_opt  = future.get();
    if (not record_ptr_opt) {

      ExecEnv::log().warn("AnalysisTranscriptSequence; problem analyzing transcript: {}", transcipt_id);
      // Skip the statistics
      continue;

    }
    // Reduce code noise.
    auto const& record_ptr = record_ptr_opt.value();
    auto const modified_view = record_ptr->modified().getView();
    auto const& genome_id = record_ptr->genome();
    std::string seq_label = TranscriptSequenceRecord::generateSequenceLabel(modified_view);

    if (transcript_map_.contains(seq_label)) {

      auto& [seq_view, sequence_record] = *transcript_map_.find(seq_label);
      if (not sequence_record.addSequenceRecord(record_ptr)) {

        ExecEnv::log().warn("Cannot add transcript record for genome: {}", genome_id);

      }

    } else {

      auto [insert_iter, result] = transcript_map_.try_emplace(seq_label, TranscriptSequenceRecord(record_ptr));
      if (not result) {

        ExecEnv::log().error("Unable to insert sequence record_ptr for genome: {}", genome_id);

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
  auto& [original_sequence, original_validity, original_amino_size] = original_opt.value();
  auto& [modified_sequence, modified_validity, modified_amino_size] = modified_opt.value();

  auto modify_record_ptr = std::make_shared<const TranscriptModifyRecord>(genome_id,
                                                                          transcript_ptr_,
                                                                          std::move(original_sequence),
                                                                          original_validity,
                                                                          std::move(modified_sequence),
                                                                          modified_validity);
  return modify_record_ptr;

}


std::string kga::TranscriptSequenceRecord::generateSequenceLabel(const DNA5SequenceCodingView& seq_view) {

  size_t sequence_tag = std::hash<std::string_view>{}(seq_view.getStringView());
  return std::format(SEQUENCE_LABEL_FORMAT_, sequence_tag);

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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

  report_stream << "Tag" << REPORT_FIELD_
                << "Genomes" << REPORT_FIELD_
                << "TranscriptCount" << REPORT_FIELD_
                << "Transcripts" << REPORT_FIELD_
                << "Text" << '\n';

  for (auto const& [seq_label, record] : transcript_map_) {

    std::string transcript_list;
    for (auto const& [transcript_id, trans_vector] : record.getTranscripts()) {

      if (not transcript_list.empty()) {

        transcript_list += REPORT_SUBFIELD_;

      }

      transcript_list += transcript_id;
      transcript_list += std::format("({})", trans_vector.size());

    }

    auto modified_sub_view = record.getModifiedView().getTruncate(SUB_VIEW_INTERVAL_);

    report_stream << seq_label << REPORT_FIELD_
                  << record.getGenomes().size() << REPORT_FIELD_
                  << record.getTranscripts().size() << REPORT_FIELD_
                  << transcript_list << REPORT_FIELD_
                  << modified_sub_view.getStringView() << '\n';

  }

}

void kga::AnalysisTranscriptSequence::mergeMap(const TranscriptMap& transcript_map) {

  for (auto const& [modified_view, merge_record] : transcript_map) {

    if (transcript_map_.contains(modified_view)) {

      auto& [map_view, map_record] = *transcript_map_.find(modified_view);
      for (auto const& [merge_genome, merge_record_ptr] : merge_record.getGenomes()) {

        if (not map_record.addSequenceRecord(merge_record_ptr)) {

          ExecEnv::log().warn("Cannot merge transcript record for genome: {}", merge_genome);

        }

      }

    } else {

      auto [insert_iter, result] = transcript_map_.try_emplace(modified_view, merge_record);
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

      auto [inter_iter, insert_result] = analysis_map_.try_emplace(transcript_id, transcript_ptr);
      if (not insert_result) {

        ExecEnv::log().warn("Unable to creat duplicate analysis for transcript: {}", transcript_id);

      }

    }

  } // For all genes.

}

void kga::AnalysisTranscriptFamily::performFamilyAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr) {


  for (auto& [transcript_id, transcript_analysis] : analysis_map_) {

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

  for (const auto& [transcript_id, transcript_analysis] : analysis_map_) {

    transcript_analysis.printReport(report_directory);

  }

  if (not analysis_map_.empty()) {

    auto merged_analysis = generateTotal();
    merged_analysis.printReport("MergedTotal", report_directory);

  }

  GenomeTranscriptAnalysis genome_analysis;
  if (not genome_analysis.createTranscriptMap(getMap())) {

    ExecEnv::log().warn("Transcript analysis index by genome failed");

  }

  genome_analysis.printGenomeReport(report_directory);

}

kga::AnalysisTranscriptSequence kga::AnalysisTranscriptFamily::generateTotal() const {

  // Copy the first view
  if (analysis_map_.empty()) {

    ExecEnv::log().critical("Unexpected empty analysis vector");

  }

  auto const& [first_transcript, first_trans_analysis] = *analysis_map_.begin();
  AnalysisTranscriptSequence merged_sequence(first_trans_analysis);
  // Vector will preserve sort order, drop first view
  for (auto const& [transcript_id, trans_sequence] : std::ranges::drop_view{analysis_map_, 1}) {

    merged_sequence.mergeMap(trans_sequence.transcriptMap());

  }

  return merged_sequence;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kga::GenomeTranscriptAnalysis::processTranscriptMap(const AnalysisTranscriptMap& transcript_map) {

  bool result{true};
  for (auto const& [transcript_id, transcript_record] : transcript_map) {

    const TranscriptMap& record_transcript_map = transcript_record.transcriptMap();
    for (auto const& [sequence_label, sequence_record] : record_transcript_map) {

      for (auto const& [genome_id, transcript_modify] : sequence_record.getGenomes()) {

        if (genome_map_.contains(genome_id)) {

          auto& [genome_key, genome_trans_map] = *genome_map_.find(genome_id);
          auto [insert_iter, insert_result] = genome_trans_map.try_emplace(transcript_id, transcript_modify);
          if (not insert_result) {

            ExecEnv::log().warn("Duplicate modify record for transcript: {}, genome: {}", transcript_id, genome_id);
            result = false;

          }

        } else {

          FeatureTranscriptMap feature_map {{transcript_id, transcript_modify}};
          auto [insert_iter, insert_result] = genome_map_.try_emplace(genome_id, feature_map);
          if (not insert_result) {

            ExecEnv::log().error("Duplicate modify record for transcript: {}, genome: {}", transcript_id, genome_id);
            result = false;

          }

        }

      }

    }

  }

  return result;

}

// Ensure all genome records contain the same modified transcripts.
std::optional<std::vector<kgl::FeatureIdent_t>> kga::GenomeTranscriptAnalysis::checkGenomeMap() const {

  if (genome_map_.empty()) {

    ExecEnv::log().warn("Transcript map contains no genomes");
    return std::nullopt;

  }

  auto const& [first_genome, first_transcript_map] = *genome_map_.begin();
  // Get the transcript vector.
  std::vector<kgl::FeatureIdent_t> transcript_vector;
  for (auto const& [feature_id, transcript_modify] : first_transcript_map) {

    transcript_vector.push_back(feature_id);

  }

  for (auto const& [genome, transcript_map] : std::ranges::drop_view{genome_map_, 1}) {

    // Check number of transcripts.
    if (transcript_map.size() == transcript_vector.size()) {

      for (auto const& [map_pair, transcript_id] : std::ranges::zip_view(transcript_map, transcript_vector)) {

        auto const& [map_transcript, map_transcript_modify] = map_pair;
        if (transcript_id != map_transcript) {

          ExecEnv::log().warn("Genome: {}, map transcript id: {} expected transcript id: {}",
                              genome, map_transcript, transcript_id);
          return std::nullopt;

        }

      }

    } else {

      ExecEnv::log().warn("Genome: {}, transcript count: {} expected transcript count: {}",
                          genome, transcript_map.size(), transcript_vector.size());
      return std::nullopt;

    }

  }

  return transcript_vector;

}

bool kga::GenomeTranscriptAnalysis::createTranscriptMap(const AnalysisTranscriptMap& transcript_map) {

  bool result{true};

  if (not processTranscriptMap( transcript_map)) {

    ExecEnv::log().warn("Process transcript map index by genome failed");
    result = false;

  }

  auto transcript_vec_opt = checkGenomeMap();
  if (not transcript_vec_opt) {

    ExecEnv::log().warn("Integrity of the genome map failed");
    result = false;

  }

  return result;

}

void kga::GenomeTranscriptAnalysis::printGenomeReport(const std::string& report_directory) const {

  // Select distance metric.
  CodingDistanceMetric dna_distance_metric{LevenshteinGlobalCoding};

  auto file_path = REPORT_PREFIX_ + REPORT_EXT_;
  file_path = kel::Utility::filePath(file_path, report_directory);

  std::ofstream report_stream(file_path);
  if (not report_stream.good()) {

    ExecEnv::log().warn("GenomeTranscriptAnalysis; Could not open file: {}", file_path);
    return;

  }

  auto transcript_vec_opt = checkGenomeMap();
  if (not transcript_vec_opt) {

    ExecEnv::log().error("Integrity of the genome map failed");
    return;

  }
  auto const& transcript_vector = transcript_vec_opt.value();

  report_stream << "Genomes";
  for (auto const& transcript_id : transcript_vector) {

    report_stream << REPORT_FIELD_
                  << transcript_id
                  << REPORT_FIELD_
                  << "Distance";

  }
  report_stream << '\n';

  for (auto const& [genome_id, transcript_map] : genome_map_) {

    report_stream << genome_id;
    for (auto const& [transcript_id, transcript_modify] : transcript_map) {

      auto modify_transcript_view = transcript_modify->modified().getView();
      auto sequence_label = TranscriptSequenceRecord::generateSequenceLabel(modify_transcript_view);
      double distance = dna_distance_metric(transcript_modify->reference(), transcript_modify->modified());
      report_stream << REPORT_FIELD_
                    << sequence_label
                    << REPORT_FIELD_
                    << distance;

    }
    report_stream << '\n';

  }

}
