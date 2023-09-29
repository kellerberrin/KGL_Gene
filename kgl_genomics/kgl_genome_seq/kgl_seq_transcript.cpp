//
// Created by kellerberrin on 16/09/23.
//

#include "kgl_seq_transcript.h"
#include "kgl_seq_offset.h"
#include "kel_interval_set.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


std::pair<kgl::SequenceStats, bool>
kgl::SequenceTranscript::createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                                const std::shared_ptr<const ContigReference>& contig_reference_ptr,
                                                const OpenRightUnsigned& sequence_interval) {

  // Return filtered variants adjusted for duplicate variants and upstream deletes.
  auto [interval_map, non_unique_count, upstream_deleted] = MutationOffset::getCanonicalVariants(contig_variant_ptr, sequence_interval);
  size_t map_size = interval_map.size() + non_unique_count + upstream_deleted;

  auto adjusted_offset_ptr = std::make_shared<AdjustedSequenceInterval>(sequence_interval);
  if (not adjusted_offset_ptr->processVariantMap(interval_map)) {

    ExecEnv::log().warn("Problem updating interval: {}, variant contig: {}, reference contig_ref_ptr: {}",
                        sequence_interval.toString(),
                        contig_variant_ptr->contigId(),
                        contig_reference_ptr->contigId());
    return {{}, false};

  }
  if (not adjusted_sequence_.updateSequence(contig_reference_ptr,
                                            sequence_interval,
                                            adjusted_offset_ptr->indelModifyMap())) {

    ExecEnv::log().warn("Problem updating sequence: {}, variant contig: {}, reference contig_ref_ptr: {}",
                        sequence_interval.toString(),
                        contig_variant_ptr->contigId(),
                        contig_reference_ptr->contigId());
    return {{}, false};

  }

  SequenceStats sequence_stats;
  sequence_stats.map_size_ = map_size;
  sequence_stats.non_unique_count_ = non_unique_count;
  sequence_stats.upstream_deleted_ = upstream_deleted;

  return { sequence_stats, true};

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::concatModifiedSequences(const std::vector<OpenRightUnsigned>& interval_vector) const {

  // Sort the intervals.
  IntervalSetLower interval_set;
  for (auto const& interval : interval_vector) {

    auto [insert_iter, result] = interval_set.insert(interval);
    if (not result) {

      ExecEnv::log().warn("Interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    auto modified_sequence_opt = adjusted_sequence_.modifiedSubSequence(sub_interval);
    if (not modified_sequence_opt) {

      ExecEnv::log().warn("Unable to generate modified sequence for interval: {}", sub_interval.toString());
      return std::nullopt;

    }

    bool result = concatenated_sequence.append(modified_sequence_opt.value());
    if (not result) {

      ExecEnv::log().warn("Unable to concatenate modified sequence for interval: {}", sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::concatOriginalSequences(const std::vector<OpenRightUnsigned>& interval_vector) const {

  // Sort the intervals.
  IntervalSetLower interval_set;
  for (auto const& interval : interval_vector) {

    auto [insert_iter, result] = interval_set.insert(interval);
    if (not result) {

      ExecEnv::log().warn("SequenceTranscript::concatOriginalSequences; interval: {} has a duplicate lower()", interval.toString());

    }

  }

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : interval_set) {

    auto modified_sequence_opt = adjusted_sequence_.originalSubSequence(sub_interval);
    if (not modified_sequence_opt) {

      ExecEnv::log().warn("Unable to generate original sequence for interval: {}", sub_interval.toString());
      return std::nullopt;

    }

    bool result = concatenated_sequence.append(modified_sequence_opt.value());
    if (not result) {

      ExecEnv::log().warn("Unable to concatenate modified sequence for interval: {}", sub_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated_sequence;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getModifiedGene(const GeneIntervalStructure& gene_interval,
                                                                                const FeatureIdent_t& transcript_id) const {

  auto find_iter = gene_interval.codingTranscripts().find(transcript_id);
  if (find_iter == gene_interval.codingTranscripts().end()) {

    ExecEnv::log().warn("Could not find transcript: {} for gene: {}", transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  auto [trans_id, exon_set] = *find_iter;
  std::vector<OpenRightUnsigned> exon_vector{exon_set.begin(), exon_set.end()};

  auto sequence_opt = concatModifiedSequences(exon_vector);
  if (not sequence_opt) {

    ExecEnv::log().warn("Could not modify transcript: {} for gene: {}", transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  return sequence_opt;

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getOriginalGene(const GeneIntervalStructure& gene_interval,
                                                                                const FeatureIdent_t& transcript_id) const {

  auto find_iter = gene_interval.codingTranscripts().find(transcript_id);
  if (find_iter == gene_interval.codingTranscripts().end()) {

    ExecEnv::log().warn("Could not find transcript: {} for gene: {}", transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  auto [trans_id, intron_set] = *find_iter;
  std::vector<OpenRightUnsigned> exon_vector{intron_set.begin(), intron_set.end()};

  auto sequence_opt = concatOriginalSequences(exon_vector);
  if (not sequence_opt) {

    ExecEnv::log().warn("Could not modify transcript: {} for gene: {}", transcript_id, gene_interval.getGene()->id());
    return std::nullopt;

  }

  return sequence_opt;

}
