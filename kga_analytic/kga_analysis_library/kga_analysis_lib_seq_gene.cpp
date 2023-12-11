//
// Created by kellerberrin on 17/11/23.
//

#include "kel_workflow_threads.h"
#include "kga_analysis_lib_seq_gene.h"
#include "kgl_sequence_motif.h"
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
                << "PfPol_III_ABox" << REPORT_FIELD_
                << "PfPol_III_BBox" << REPORT_FIELD_
                << "Sequence" << '\n';

  for (auto const& [seq_label, record] : transcript_map_) {

    std::string transcript_list;
    for (auto const& [transcript_id, trans_vector] : record.getTranscripts()) {

      if (not transcript_list.empty()) {

        transcript_list += REPORT_SUBFIELD_;

      }

      transcript_list += transcript_id;
      transcript_list += std::format("({})", trans_vector.size());

    }

    auto modified_sub_view = record.getModifiedView().getIntersection(SUB_VIEW_INTERVAL_);
    report_stream << seq_label << REPORT_FIELD_
                  << record.getGenomes().size() << REPORT_FIELD_
                  << record.getTranscripts().size() << REPORT_FIELD_
                  << transcript_list << REPORT_FIELD_
                  << SearchSequence::PfPolymerase_III_ABox(record.getModifiedView()).size() << REPORT_FIELD_
                  << SearchSequence::PfPolymerase_III_BBox(record.getModifiedView()).size() << REPORT_FIELD_
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


