//
// Created by kellerberrin on 9/12/23.
//

#include "kga_analysis_lib_seq_stats.h"
#include "kgl_sequence_distance_impl.h"
#include "kgl_distance_tree_upgma.h"
#include "kgl_sequence_node_view.h"
#include "kgl_classification_tree.h"


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


void kga::AnalysisTranscriptFamily::printAllReports(const std::string& analysis_directory,
                                                    const std::string& analysis_sub_directory,
                                                    const std::map<GenomeId_t, std::string>& annotation_map) const {

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

  genome_analysis.printGenomeReport(report_directory, annotation_map);


}


void  kga::AnalysisTranscriptFamily::createClassificationTree(const std::string& analysis_directory,
                                                              const std::string& analysis_sub_directory,
                                                              const std::map<GenomeId_t, std::string>& annotations,
                                                              const std::string& selection_text) const {

  auto report_directory = Utility::appendPath(analysis_sub_directory, analysis_directory);

  // Recreate the report subdirectory.
  if (not Utility::directoryRenew(report_directory)) {

    ExecEnv::log().warn("AnalysisTranscriptFamily; cannot re-create report directory: {}", report_directory);

  }

  GenomeTranscriptAnalysis genome_analysis;
  if (not genome_analysis.createTranscriptMap(analysis_map_)) {

    ExecEnv::log().warn("Transcript analysis index by genome failed");

  }

  for (auto const& [transcript_id, transcript_analysis] : analysis_map_) {


    genome_analysis.createClassificationTree( report_directory, transcript_id, annotations, selection_text);


  }

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

void kga::GenomeTranscriptAnalysis::printGenomeReport(const std::string& report_directory,
                                                      const std::map<GenomeId_t, std::string>& annotations) const {

  // Select parentDistance metric.
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

  report_stream << "Genomes" << REPORT_FIELD_ << "Annotation";
  for (auto const& transcript_id : transcript_vector) {

    report_stream << REPORT_FIELD_
                  << transcript_id
                  << REPORT_FIELD_
                  << "Distance";

  }
  report_stream << REPORT_FIELD_ << "TotalDistance";
  report_stream << '\n';

  for (auto const& [genome_id, transcript_map] : genome_map_) {

    report_stream << genome_id << REPORT_FIELD_;
    if (annotations.contains(genome_id)) {

      auto const& [genome_key, genome_annotation] = *annotations.find(genome_id);
      report_stream << genome_annotation;

    } else {

      report_stream << "n/a";

    }

    size_t total_distance{0};
    for (auto const& [transcript_id, transcript_modify] : transcript_map) {

      auto modify_transcript_view = transcript_modify->modified().getView();
      auto sequence_label = TranscriptSequenceRecord::generateSequenceLabel(modify_transcript_view);
      double distance = dna_distance_metric(transcript_modify->reference(), transcript_modify->modified());
      total_distance += distance;
      report_stream << REPORT_FIELD_
                    << sequence_label
                    << REPORT_FIELD_
                    << distance;

    }
    report_stream << REPORT_FIELD_ << total_distance;
    report_stream << '\n';

  }

}

void kga::GenomeTranscriptAnalysis::createClassificationTree( const std::string& newick_directory,
                                                              const std::string& transcript_id,
                                                              const std::map<GenomeId_t, std::string>& annotations,
                                                              const std::string& selection_text) const {



  TreeNodeVector tree_node_vector;
  CodingDistanceMetricView coding_distance_metric_view{LevenshteinLocalCodingView};
  for (auto const& [genome_id, transcript_map] : genome_map_) {

    if (transcript_map.contains(transcript_id)) {

      std::string node_tag = genome_id;
      if (annotations.contains(genome_id)) {

        auto const& [genome_key, genome_annotation] = *annotations.find(genome_id);
        node_tag += "_" + genome_annotation;

      }

      if (not selection_text.empty() and not node_tag.contains(selection_text)) {

        // Skip this record.
        continue;

      }

      // Note that the modified sequences are held in 'genome_map_' and will exist in memory for the tree analysis.
      auto const& [transcript_key, transcript_record] = *transcript_map.find(transcript_id);
      // For performance reasons, sequence view is analysed in the tree, therefore the underlying sequence must be in memory during analysis.
      auto const& modified_transcript_view = transcript_record->modified().getView();

      auto coding_sequence_ptr = std::make_shared<CodingSequenceViewNode>(modified_transcript_view, node_tag, coding_distance_metric_view);
      tree_node_vector.push_back(coding_sequence_ptr);

    } // Transcript

  }

  ExecEnv::log().info("Tree nodes: {} created for: {}", tree_node_vector.size(), transcript_id);

  DistanceTreeUPGMA upgma_distance;
  // Add the parentDistance vector.
  upgma_distance.addDistanceMap(tree_node_vector);
  // Calculate.
  auto root_node_vector = upgma_distance.calculateTree();
  // Report Results.
  ClassificationTree upgma_tree(root_node_vector);

  auto file_path = transcript_id + "_" + NEWICK_SUFFIX_ + NEWICK_EXT_;
  file_path = kel::Utility::filePath(file_path, newick_directory);

  if (not upgma_tree.writeNewick(file_path)) {

    ExecEnv::log().warn("createClassificationTree; Could not open file: {}", file_path);

  }

}
