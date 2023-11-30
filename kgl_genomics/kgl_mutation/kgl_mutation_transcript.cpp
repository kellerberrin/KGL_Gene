//
// Created by kellerberrin on 16/09/23.
//

#include "kgl_mutation_transcript.h"
#include "kgl_mutation_variant_filter.h"
#include "kgl_sequence_codon.h"


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


std::pair<kgl::FilteredVariantStats, bool>
kgl::SequenceTranscript::createModifiedSequence(const std::shared_ptr<const ContigDB>& contig_variant_ptr,
                                                SeqVariantFilterType filter_type) {


  // Request (not guaranteed) the extension of the transcript interval by a 5_prime buffer and a 3_prime_buffer.
  prime_5_extend_ = transcript_ptr_->prime5Region(PRIME_5_BUFFER_); // Store the actual extension.
  prime_3_extend_ = transcript_ptr_->prime3Region(PRIME_3_BUFFER_);
  auto extended_transcript = transcript_ptr_->extendInterval(PRIME_5_BUFFER_, PRIME_3_BUFFER_);
  // Return filtered variants adjusted for duplicate variants and upstream deletes.
  SequenceVariantFilter filtered_variants(contig_variant_ptr, extended_transcript, filter_type);

  const auto& contig_reference_ptr = transcript_ptr_->contig();
  if (not adjusted_sequence_.updateSequence(contig_reference_ptr, filtered_variants)) {

    ExecEnv::log().warn("Problem updating sequence: {}, variant contig: {}, reference contig_ref_ptr: {}",
                        extended_transcript.toString(),
                        contig_variant_ptr->contigId(),
                        contig_reference_ptr->contigId());
    return {{}, false};

  }

  return { filtered_variants.filterStatistics(), true};

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


std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>> kgl::SequenceTranscript::getModifiedValidity() const {

  auto modified_linear_opt = getModifiedLinear();
  if (not modified_linear_opt) {

    return std::nullopt;

  }

  return getValidity(std::move(modified_linear_opt.value()));

}


std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>> kgl::SequenceTranscript::getModifiedAdjustedValidity() const {

  auto modified_linear_opt = getModifiedAdjusted();
  if (not modified_linear_opt) {

    return std::nullopt;

  }

  return getValidity(std::move(modified_linear_opt.value()));

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getOriginalLinear() const {

  // Extract the modified sequences and concatenate them.
  DNA5SequenceLinear concatenated_sequence;
  for (auto const& sub_interval : transcript_ptr_->getExonIntervals()) {

    auto modified_sequence_opt = adjusted_sequence_.originalSubSequence(sub_interval);
    if (not modified_sequence_opt) {

      ExecEnv::log().warn("Unable to generate reference sequence for interval: {}", sub_interval.toString());
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


std::optional<std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>> kgl::SequenceTranscript::getOriginalValidity() const {

  auto original_linear_opt = getOriginalLinear();
  if (not original_linear_opt) {

    return std::nullopt;

  }

  return getValidity(std::move(original_linear_opt.value()));

}


std::optional<kgl::DNA5SequenceLinear> kgl::SequenceTranscript::getModifiedAdjusted() const {


  auto modified_linear_opt = getModifiedLinear();
  if (not modified_linear_opt) {

    return std::nullopt;

  }
  auto modified_linear(std::move(modified_linear_opt.value()));
  size_t original_length = modified_linear.length();

  // Check mod3
  size_t mod3_remainder = Codon::codonRemainder(original_length);
  size_t adjusted_size = Codon::CODON_SIZE - mod3_remainder;
  size_t extend_size{0};
  bool adjust_flag{false};

  if (mod3_remainder != 0) {

    auto prime3_sequence_opt = adjusted_sequence_.modifiedSubSequence(prime_3_extend_);
    if (not prime3_sequence_opt) {

      ExecEnv::log().warn("Unable to generate extended modified sequence for interval: {}", prime_3_extend_.toString());
      return std::nullopt;

    }
    const auto& prime3_sequence = prime3_sequence_opt.value();

    if (adjusted_size <= prime3_sequence.length()) {

      OpenRightUnsigned extend_interval{0, adjusted_size};
      auto extend_sequence_opt = prime3_sequence.subSequence(extend_interval);
      if (not extend_sequence_opt) {

        ExecEnv::log().warn("Unable to extract extend sequence for interval: {}", extend_interval.toString());
        return std::nullopt;

      }

      bool result = modified_linear.append(extend_sequence_opt.value());
      if (not result) {

        ExecEnv::log().warn("Unable to concatenate modified sequence for interval: {}", extend_interval.toString());
        return std::nullopt;

      }

      adjust_flag = true;
      extend_size = extend_sequence_opt.value().length();

    }

  }

  if (Codon::codonRemainder(modified_linear.length()) != 0 and adjust_flag) {

    ExecEnv::log().warn("Adjusted sequence not mod3: original length: {}, reference mod3: {}, adjusted size: {}, extend size: {}, modified length: {}, mod3 remainder adjusted: {}",
                        original_length,
                        mod3_remainder,
                        adjusted_size,
                        extend_size,
                        modified_linear.length(),
                        Codon::codonRemainder(modified_linear.length()));

  }

  return modified_linear;

}


std::tuple<kgl::DNA5SequenceCoding, kgl::CodingSequenceValidity, size_t>
kgl::SequenceTranscript::getValidity(DNA5SequenceLinear&& linear_coding) const {

  auto modified_coding = linear_coding.codingSequence(transcript_ptr_->strand());

  CodingSequenceValidity sequence_validity;
  size_t amino_size{0};
  if (transcript_ptr_->codingType() == TranscriptionSequenceType::PROTEIN) {

    const auto& contig_ref_ptr = transcript_ptr_->getGene()->contig_ref_ptr();
    auto [validity, sequence_size] = contig_ref_ptr->codingProteinSequenceSize(modified_coding);
    sequence_validity = validity;
    amino_size = sequence_size;

  } else {

    sequence_validity = CodingSequenceValidity::NCRNA;

  }

  return {std::move(modified_coding), sequence_validity, amino_size};

}



