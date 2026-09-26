//
// Consensus library (glm-refactor polish of deepseek-refactor).
//
// Module 4 implementation: gene/transcript client.
//

#include "kgl_consensus_transcript.h"
#include "kgl_sequence_codon.h"
#include "kel_exec_env.h"

#include <algorithm>

namespace kgl = kellerberrin::genome;

kgl::SequenceTranscript::SequenceTranscript(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                            std::shared_ptr<const TranscriptionSequence> transcript_ptr,
                                            SeqVariantFilterType filter_type)
    : transcript_ptr_(std::move(transcript_ptr)) {

  createModifiedSequence(contig_variant_ptr, filter_type);

}

void kgl::SequenceTranscript::createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                                     SeqVariantFilterType filter_type) {

  sequence_status_ = false;

  if (not transcript_ptr_ or not contig_variant_ptr) {

    ExecEnv::log().warn("SequenceTranscript; null transcript or variant contig");
    return;

  }

  // The extended transcript window: 5' buffer (currently 0) and 3' mod3 buffer (200).
  prime_5_extend_ = transcript_ptr_->prime5Region(PRIME_5_BUFFER_);
  prime_3_extend_ = transcript_ptr_->prime3Region(PRIME_3_BUFFER_);
  window_ = transcript_ptr_->extendInterval(PRIME_5_BUFFER_, PRIME_3_BUFFER_);

  // Module 1: window-scoped selection. The statistics are per extended-transcript
  // window (the client contract); the returned set is policy-filtered, resolved to
  // one variant per locus and shadow-pruned.
  auto selected_opt = selectWindowVariants(contig_variant_ptr, window_, filter_type);
  if (not selected_opt) {

    ExecEnv::log().warn("SequenceTranscript; variant selection failed: {}", toString(selected_opt.error()));
    return;

  }
  // Module 2: global edit schedule and offset accounting.
  // The resolved schedule owns the upstream_deleted_ statistic: it is taken
  // after module 1's shadow pruning plus module 2's delete-union pruning, so
  // edits erased by merged/adjacent delete spans are also counted (single-sourced).
  // In the production path module 1's ContigUpstreamFilter already removed
  // overlapping-delete conflicts, so the module 2 increment mostly fires for
  // adjacent deletes and for directly constructed SelectedVariants (tests, probes).
  auto applied_opt = resolveVariants(*selected_opt);
  if (not applied_opt) {

    ExecEnv::log().warn("SequenceTranscript; variant resolution failed: {}", toString(applied_opt.error()));
    return;

  }
  selected_opt->stats.upstream_deleted_ = applied_opt->stats.upstream_deleted_;
  filter_stats_ = selected_opt->stats;

  auto accounting_opt = OffsetAccounting::build(*applied_opt);
  if (not accounting_opt) {

    ExecEnv::log().warn("SequenceTranscript; offset accounting failed: {}", toString(accounting_opt.error()));
    return;

  }

  // Module 3: materialise the consensus over the extended window.
  const auto& contig_ref_ptr = transcript_ptr_->contig();
  auto reference_opt = contig_ref_ptr->sequence().subSequence(window_);
  if (not reference_opt) {

    ExecEnv::log().warn("SequenceTranscript; cannot extract reference window: {}", window_.toString());
    return;

  }
  original_window_ = std::move(reference_opt.value());

  auto consensus_opt = buildConsensus(original_window_, window_, *accounting_opt);
  if (not consensus_opt) {

    ExecEnv::log().warn("SequenceTranscript; consensus build failed: {}", toString(consensus_opt.error()));
    return;

  }
  consensus_.emplace(std::move(consensus_opt.value()));
  sequence_status_ = true;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::modifiedLinear() const {

  if (not consensus_) {

    return std::nullopt;

  }

  DNA5SequenceLinear concatenated;
  for (auto const& exon_interval : transcript_ptr_->getExonIntervals()) {

    // Exons outside the extended window cannot contribute (the window covers the gene).
    auto slice_opt = consensus_->slice(exon_interval);
    if (not slice_opt) {

      ExecEnv::log().warn("SequenceTranscript; exon interval {} outside extended window {}",
                          exon_interval.toString(), window_.toString());
      return std::nullopt;

    }
    if (not concatenated.append(slice_opt.value())) {

      ExecEnv::log().warn("SequenceTranscript; cannot concatenate exon interval: {}", exon_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::originalLinear() const {

  DNA5SequenceLinear concatenated;
  for (auto const& exon_interval : transcript_ptr_->getExonIntervals()) {

    if (not window_.containsInterval(exon_interval)) {

      ExecEnv::log().warn("SequenceTranscript; exon interval {} outside window {}",
                          exon_interval.toString(), window_.toString());
      return std::nullopt;

    }
    const OpenRightUnsigned local{exon_interval.lower() - window_.lower(),
                                  exon_interval.upper() - window_.lower()};
    auto sub_opt = original_window_.subSequence(local);
    if (not sub_opt) {

      ExecEnv::log().warn("SequenceTranscript; cannot extract original exon interval: {}", exon_interval.toString());
      return std::nullopt;

    }
    if (not concatenated.append(sub_opt.value())) {

      ExecEnv::log().warn("SequenceTranscript; cannot concatenate original exon interval: {}", exon_interval.toString());
      return std::nullopt;

    }

  }

  return concatenated;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getModifiedLinear() const {

  return sequence_status_ ? modifiedLinear() : std::nullopt;

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getOriginalLinear() const {

  return sequence_status_ ? originalLinear() : std::nullopt;

}

std::optional<kgl::DNA5SequenceCoding> kgl::SequenceTranscript::getModifiedCoding() const {

  auto linear_opt = getModifiedLinear();
  if (not linear_opt) {

    return std::nullopt;

  }
  return linear_opt.value().codingSequence(transcript_ptr_->strand());

}

std::optional<kgl::DNA5SequenceCoding> kgl::SequenceTranscript::getOriginalCoding() const {

  auto linear_opt = getOriginalLinear();
  if (not linear_opt) {

    return std::nullopt;

  }
  return linear_opt.value().codingSequence(transcript_ptr_->strand());

}

std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>>
kgl::SequenceTranscript::getModifiedValidity() const {

  auto linear_opt = getModifiedLinear();
  if (not linear_opt) {

    return std::nullopt;

  }
  return getValidity(std::move(linear_opt.value()));

}

std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>>
kgl::SequenceTranscript::getModifiedAdjustedValidity() const {

  auto linear_opt = getModifiedAdjusted();
  if (not linear_opt) {

    return std::nullopt;

  }
  return getValidity(std::move(linear_opt.value()));

}

std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>>
kgl::SequenceTranscript::getOriginalValidity() const {

  auto linear_opt = getOriginalLinear();
  if (not linear_opt) {

    return std::nullopt;

  }
  return getValidity(std::move(linear_opt.value()));

}

std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getModifiedAdjusted() const {

  auto linear_opt = modifiedLinear();
  if (not linear_opt) {

    return std::nullopt;

  }
  auto modified_linear = std::move(linear_opt.value());
  const size_t original_length = modified_linear.length();
  const size_t remainder = Codon::codonRemainder(original_length);
  const size_t adjusted_size = (Codon::CODON_SIZE - remainder) % Codon::CODON_SIZE;

  if (remainder != 0 and consensus_) {

    auto prime3_opt = consensus_->slice(prime_3_extend_);
    if (not prime3_opt) {

      ExecEnv::log().warn("SequenceTranscript; cannot extract 3' buffer: {}", prime_3_extend_.toString());
      return std::nullopt;

    }
    const auto& prime3 = prime3_opt.value();
    if (adjusted_size <= prime3.length()) {

      auto extend_opt = prime3.subSequence({0, adjusted_size});
      if (not extend_opt or not modified_linear.append(extend_opt.value())) {

        ExecEnv::log().warn("SequenceTranscript; cannot append 3' mod3 bases");
        return std::nullopt;

      }

    }

  }

  if (Codon::codonRemainder(modified_linear.length()) != 0 and remainder != 0) {

    ExecEnv::log().warn("SequenceTranscript; adjusted sequence not mod3: original: {}, adjusted: {}",
                        original_length, modified_linear.length());

  }

  return modified_linear;

}

std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>
kgl::SequenceTranscript::getValidity(DNA5SequenceLinear&& linear_coding) const {

  auto modified_coding = linear_coding.codingSequence(transcript_ptr_->strand());
  CodingSequenceValidity validity{CodingSequenceValidity::NCRNA};
  size_t amino_size{0};

  if (transcript_ptr_->codingType() == TranscriptionSequenceType::PROTEIN) {

    const auto& contig_ref_ptr = transcript_ptr_->getGene()->contig_ref_ptr();
    auto [protein_validity, sequence_size] = contig_ref_ptr->codingProteinSequenceSize(modified_coding);
    validity = protein_validity;
    amino_size = sequence_size;

  }

  return {std::move(modified_coding), validity, amino_size};

}
