//
// Created by kellerberrin on 16/09/23.
//

#include "kgl_mutation_aggregation.h"
#include "kel_interval_set.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


std::optional<kgl::DNA5SequenceLinear>
    kgl::SequenceAggregation::concatModifiedSequences( const AdjustedSequence& adjusted_sequence,
                                                       const std::vector<OpenRightUnsigned>& interval_vector) {

  // Sort the intervals.
  IntervalSetLower interval_set;
  for (auto const& interval : interval_vector) {

    auto [insert_iter, result] = interval_set.insert(interval);
    if (not result) {

      ExecEnv::log().warn("SequenceAggregation::concatModifiedSequences; interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    auto modified_sequence_opt = adjusted_sequence.modifiedSubSequence(sub_interval);
    if (not modified_sequence_opt) {

      ExecEnv::log().warn("SequenceAggregation::concatModifiedSequences; unable to generate modified sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

    bool result = concatenated_sequence.append(modified_sequence_opt.value());
    if (not result) {

      ExecEnv::log().warn("SequenceAggregation::concatModifiedSequences; unable to concatenate modified sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}


std::optional<kgl::DNA5SequenceLinear>
kgl::SequenceAggregation::concatOriginalSequences( const AdjustedSequence& adjusted_sequence,
                                                   const std::vector<OpenRightUnsigned>& interval_vector) {

  // Sort the intervals.
  IntervalSetLower interval_set;
  for (auto const& interval : interval_vector) {

    auto [insert_iter, result] = interval_set.insert(interval);
    if (not result) {

      ExecEnv::log().warn("SequenceAggregation::concatOriginalSequences; interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    auto modified_sequence_opt = adjusted_sequence.originalSubSequence(sub_interval);
    if (not modified_sequence_opt) {

      ExecEnv::log().warn("SequenceAggregation::concatOriginalSequences; unable to generate original sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

    bool result = concatenated_sequence.append(modified_sequence_opt.value());
    if (not result) {

      ExecEnv::log().warn("SequenceAggregation::concatOriginalSequences; unable to concatenate modified sequence for interval: {}",
                          sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceAggregation::getModifiedGene( const AdjustedSequence& adjusted_sequence,
                                                                                  const GeneIntervalStructure& gene_interval,
                                                                                  const FeatureIdent_t& transcript_id) {

  auto find_iter = gene_interval.codingTranscripts().find(transcript_id);
  if (find_iter == gene_interval.codingTranscripts().end()) {

    ExecEnv::log().warn("SequenceAggregation::getModifiedGene; could not find transcript: {} for gene: {}",
                        transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  auto [trans_id, exon_set] = *find_iter;
  std::vector<OpenRightUnsigned> exon_vector{exon_set.begin(), exon_set.end()};

  auto sequence_opt = concatModifiedSequences(adjusted_sequence, exon_vector);
  if (not sequence_opt) {

    ExecEnv::log().warn("SequenceAggregation::getModifiedGene; could not modify transcript: {} for gene: {}",
                        transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  return sequence_opt;

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceAggregation::getOriginalGene( const AdjustedSequence& adjusted_sequence,
                                                                                  const GeneIntervalStructure& gene_interval,
                                                                                  const FeatureIdent_t& transcript_id) {

  auto find_iter = gene_interval.codingTranscripts().find(transcript_id);
  if (find_iter == gene_interval.codingTranscripts().end()) {

    ExecEnv::log().warn("SequenceAggregation::getOriginalGene; could not find transcript: {} for gene: {}",
                        transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  auto [trans_id, intron_set] = *find_iter;
  std::vector<OpenRightUnsigned> intron_vector{intron_set.begin(), intron_set.end()};

  auto sequence_opt = concatOriginalSequences(adjusted_sequence, intron_vector);
  if (not sequence_opt) {

    ExecEnv::log().warn("SequenceAggregation::getOriginalGene; could not modify transcript: {} for gene: {}",
                        transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  return sequence_opt;

}
