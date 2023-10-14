//
// Created by kellerberrin on 16/09/23.
//

#include "kgl_mutation_transcript.h"
#include "kgl_mutation_variant_filter.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


std::pair<kgl::SequenceStats, bool>
kgl::SequenceTranscript::createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                                SeqVariantFilterType filter_type) {

  const auto& contig_reference_ptr = transcript_ptr_->contig();

  // Return filtered variants adjusted for duplicate variants and upstream deletes.
  SequenceVariantFilter filtered_variants(contig_variant_ptr, transcript_ptr_->interval(), filter_type);
  size_t map_size = filtered_variants.preFilterVariants();

  if (not adjusted_sequence_.updateSequence(contig_reference_ptr, filtered_variants)) {

    ExecEnv::log().warn("Problem updating sequence: {}, variant contig: {}, reference contig_ref_ptr: {}",
                        transcript_ptr_->interval().toString(),
                        contig_variant_ptr->contigId(),
                        contig_reference_ptr->contigId());
    return {{}, false};

  }

  SequenceStats sequence_stats;
  sequence_stats.map_size_ = map_size;
  sequence_stats.non_unique_count_ = filtered_variants.duplicateVariants();
  sequence_stats.upstream_deleted_ = filtered_variants.downstreamDelete();

  return { sequence_stats, true};

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getModifiedLinear() const {

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : transcript_ptr_->getExonIntervals()) {

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


std::optional<kgl::DNA5SequenceCoding> kgl::SequenceTranscript::getModifiedCoding() const {

  auto modified_linear_opt = getModifiedLinear();
  if (not modified_linear_opt) {

    return std::nullopt;

  }
  const auto& modified_linear = modified_linear_opt.value();

  return modified_linear.codingSequence(transcript_ptr_->strand());

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getOriginalLinear() const {

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : transcript_ptr_->getExonIntervals()) {

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


std::optional<kgl::DNA5SequenceCoding> kgl::SequenceTranscript::getOriginalCoding() const {

  auto original_linear_opt = getOriginalLinear();
  if (not original_linear_opt) {

    return std::nullopt;

  }
  const auto& original_linear = original_linear_opt.value();

  return original_linear.codingSequence(transcript_ptr_->strand());

}
