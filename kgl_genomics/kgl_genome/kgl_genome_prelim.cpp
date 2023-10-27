//
// Created by kellerberrin on 7/11/17.
//

#include "kel_exec_env.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"

#include <ranges>


namespace kgl = kellerberrin::genome;
namespace kel = kellerberrin;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Sequence - Members
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Returns the strand adjusted feature begin (-ve is end-1) to the
// target strand adjusted feature begin (-ve is end-1)
// and returns the relative begin transcription distance as a +ve offset.
kgl::ContigOffset_t kgl::FeatureSequence::distance(const FeatureSequence& compare_feature) const {

  ContigOffset_t target_begin_offset;
  if (compare_feature.strand() == StrandSense::FORWARD) {

    target_begin_offset = compare_feature.begin();

  } else {

    target_begin_offset = compare_feature.end() - 1;

  }

  ContigOffset_t begin_offset;
  if (strand() == StrandSense::FORWARD) {

    begin_offset = begin();

  } else {

    begin_offset = end() - 1;

  }

  if (begin_offset <= target_begin_offset) {

    return target_begin_offset - begin_offset;

  }

  return begin_offset - target_begin_offset;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranscriptionSequence - Members
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kel::IntervalSetLower kgl::TranscriptionSequence::getExonIntervals() const {

  IntervalSetLower exon_intervals;
  for (const auto& [offset, feature_ptr] : transcription_feature_map_) {

    exon_intervals.insert(feature_ptr->sequence().interval());

  }

  return exon_intervals;

}

// Empty for a 1 exon gene.
kel::IntervalSetLower kgl::TranscriptionSequence::getIntronIntervals() const {

  IntervalSetLower intron_set;
  const auto exon_interval_set = getExonIntervals();

  // Iterate across a sliding window of high and low exons.
  for (auto const& [low_exon, high_exon] : std::ranges::views::pairwise(exon_interval_set)) {

    if (low_exon.disjoint(high_exon) and not low_exon.adjacent(high_exon)) {

      OpenRightUnsigned intron_interval{low_exon.upper(), high_exon.lower()};
      auto [insert_iter, result] = intron_set.insert(intron_interval);
      if (not result) {

        ExecEnv::log().warn("Unable to enter insert intron interval: {} (duplicate)", intron_interval.toString());

      }

    }

  }

  return intron_set;

}


kgl::StrandSense kgl::TranscriptionSequence::strand() const {

  return getGene()->sequence().strand();

}


// Offset of the start of the sequence - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::TranscriptionSequence::start() const {

  // Safety first.
  if (transcription_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().begin();
  }

  return transcription_feature_map_.begin()->second->sequence().begin();

}

// Offset of the end of the sequence (last nucleotide + 1) - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::TranscriptionSequence::end() const {

  // Safety first.
  if (transcription_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().end();
  }

  return transcription_feature_map_.rbegin()->second->sequence().end();

}


kel::OpenRightUnsigned kgl::TranscriptionSequence::prime5Region(ContigSize_t requested_size) const {

  // Safety first.
  if (transcription_feature_map_.empty()) {

    ExecEnv::log().error("Prime_5 sequence for gene id: {} is empty", getGene()->id());
    return {0, 0};

  }

  ContigOffset_t end_offset{0};
  ContigOffset_t begin_offset{0};
  switch(strand()) {

    case StrandSense::FORWARD:
      end_offset = transcription_feature_map_.begin()->second->sequence().begin();
      begin_offset = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(end_offset) - static_cast<SignedOffset_t>(requested_size));
      break;

    case StrandSense::REVERSE:
      begin_offset = transcription_feature_map_.rbegin()->second->sequence().end();
      end_offset = begin_offset + requested_size;
      break;

  }

  OpenRightUnsigned prime_5_interval{ begin_offset, end_offset};
  // Ensure the interval is within the contig.
  prime_5_interval = prime_5_interval.intersection(contig()->sequence().interval());

  return prime_5_interval;

}


kel::OpenRightUnsigned kgl::TranscriptionSequence::prime3Region(ContigSize_t requested_size) const {

  // Safety first.
  if (transcription_feature_map_.empty()) {

    ExecEnv::log().error("Prime_3 sequence for gene id: {} is empty", getGene()->id());
    return {0, 0};

  }

  ContigOffset_t end_offset{0};
  ContigOffset_t begin_offset{0};
  switch(strand()) {

    case StrandSense::REVERSE:
      end_offset = transcription_feature_map_.begin()->second->sequence().begin();
      begin_offset = std::max<SignedOffset_t>(0, static_cast<SignedOffset_t>(end_offset) - static_cast<SignedOffset_t>(requested_size));
      break;

    case StrandSense::FORWARD:
      begin_offset = transcription_feature_map_.rbegin()->second->sequence().end();
      end_offset = begin_offset + requested_size;
      break;

  }

  OpenRightUnsigned prime_3_interval{ begin_offset, end_offset};
  // Ensure the interval is within the contig.
  prime_3_interval = prime_3_interval.intersection(contig()->sequence().interval());

  return { prime_3_interval};

}

kel::OpenRightUnsigned kgl::TranscriptionSequence::extendInterval(ContigSize_t request_5_extend, ContigSize_t request_3_extend) const {

  // Extend the transcript by the 3_prime buffer.
  auto const transcript_interval =  interval();
  auto extended_interval = transcript_interval;

  if (request_5_extend > 0) {

    auto const prime_5_interval = prime3Region(request_5_extend);
    extended_interval = prime_5_interval.merge(extended_interval);
    if (extended_interval.empty() or extended_interval.size() <= transcript_interval.size()) {

      ExecEnv::log().warn("5 prime interval: {}, does not extend transcript interval: {}, transcript: {}",
                          prime_5_interval.toString(), transcript_interval.toString(), getParent()->id());
      return {0, 0};

    }

  }

  if (request_3_extend > 0) {

    auto prime_3_interval = prime3Region(request_3_extend);
    extended_interval = extended_interval.merge(prime_3_interval);
    if (extended_interval.empty() or extended_interval.size() <= transcript_interval.size()) {

      ExecEnv::log().warn("3 prime interval: {}, does not extend transcript interval: {}, transcript: {}",
                          prime_3_interval.toString(), transcript_interval.toString(),  getParent()->id());
      return {0, 0};

    }

  }

  return extended_interval;

}


std::shared_ptr<const kgl::ContigReference> kgl::TranscriptionSequence::contig() const {

  return getGene()->contig_ref_ptr();

}


kgl::ContigSize_t kgl::TranscriptionSequence::codingNucleotides() const {

  ContigSize_t coding_size = 0;
  // Loop through and test membership of each cds. Reminder; testing for [begin, end)
  for (const auto& cds : transcription_feature_map_) {

    coding_size += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  return coding_size;

}



void kgl::TranscriptionSequenceArray::printSequence(std::shared_ptr<const TranscriptionSequenceArray> coding_seq_ptr) {

  long vector_count = 0;
  for (const auto& sequence : coding_seq_ptr->getMap()) {

    ExecEnv::log().info("Gene: {}, begin: {}, end: {} strand: {}",
                        sequence.second->getGene()->id(),
                        sequence.second->getGene()->sequence().begin(),
                        sequence.second->getGene()->sequence().end(),
                        static_cast<char>(sequence.second->getGene()->sequence().strand()));

    ExecEnv::log().info("Parent Feature: {}, begin: {}, end: {} strand: {}",
                        sequence.second->getParent()->id(),
                        sequence.second->getParent()->sequence().begin(),
                        sequence.second->getParent()->sequence().end(),
                        static_cast<char>(sequence.second->getParent()->sequence().strand()));

    ++vector_count;

    ExecEnv::log().info("++++++++++++++ CDS Vector : {} ********************", vector_count);

    for (const auto& cds : sequence.second->getFeatureMap()) {

      ExecEnv::log().info("CDS: {}, Type: {}, begin: {}, end: {} strand: {}",
                          cds.second->id(),
                          cds.second->type(),
                          cds.second->sequence().begin(),
                          cds.second->sequence().end(),
                          static_cast<char>(cds.second->sequence().strand()));

    }

  }

}

// Assumes that all coding sequences are the same type.
kgl::TranscriptionSequenceType kgl::TranscriptionSequence::codingType() const {

  if (transcription_feature_map_.empty()) {

    return TranscriptionSequenceType::EMPTY;

  }

  auto const& [offset, feature_ptr]  = *transcription_feature_map_.begin();

  if (feature_ptr->superType() == Feature::CDS_TYPE_) {

    return TranscriptionSequenceType::PROTEIN;

  } else {

    return TranscriptionSequenceType::NCRNA;

  }

}


kgl::CodingSequenceValidity
kgl::TranscriptionSequence::checkSequenceStatus(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr) {

  if (transcript_ptr->getFeatureMap().empty()) {

    return CodingSequenceValidity::EMPTY;

  }

  if (transcript_ptr->codingType() == TranscriptionSequenceType::NCRNA) {

    return CodingSequenceValidity::NCRNA;

  }

  auto contig_ref_ptr = transcript_ptr->getGene()->contig_ref_ptr();
  auto coding_sequence_opt = contig_ref_ptr->codingSequence(transcript_ptr);
  if (not coding_sequence_opt) {

    ExecEnv::log().info("Cannot generate valid coding sequence for Gene: {}, Transcript: {}",
                        transcript_ptr->getGene()->id(),
                        transcript_ptr->getParent()->id());

    return CodingSequenceValidity::EMPTY;

  }

  DNA5SequenceCoding &coding_sequence = coding_sequence_opt.value();

  return contig_ref_ptr->checkValidCodingSequence(coding_sequence);

}


bool kgl::TranscriptionSequenceArray::insertSequence(std::shared_ptr<const TranscriptionSequence> coding_sequence_ptr) {

#ifdef CODING_SEQUENCE_ISMAP // Using a map

  auto insert = transcription_sequence_map_.insert(std::make_pair(coding_sequence_ptr->getCDSParent()->id(),
                                                           coding_sequence_ptr));

  if (not insert.second) {

    ExecEnv::log().warn("Duplicate CDS parent: {} at contig offset: {}",
                        coding_sequence_ptr->getCDSParent()->id(),
                        coding_sequence_ptr->getCDSParent()->sequence().begin());
  }

  return insert.second;

#else // Using a multimap. Same parent features permitted 

  transcription_sequence_map_.insert(std::make_pair(coding_sequence_ptr->getParent()->id(), coding_sequence_ptr));

  return true;

#endif

}

// Assumes that all coding sequences are the same type.
kgl::TranscriptionSequenceType kgl::TranscriptionSequenceArray::codingType() const {

  if (transcription_sequence_map_.empty()) {

    return TranscriptionSequenceType::EMPTY;

  }

  auto const& first_sequence = getFirst();

  return first_sequence->codingType();

}


std::shared_ptr<const kgl::TranscriptionSequence> kgl::TranscriptionSequenceArray::getFirst() const {

  if (empty()) {

    ExecEnv::log().critical("TranscriptionSequenceArray::getFirst() tried to obtain a TranscriptionSequence from empty container");

  }

  return getMap().begin()->second;

}

std::unique_ptr<const kgl::TranscriptionSequenceArray> kgl::TranscriptionSequenceArray::validTranscriptionArray() const {

  auto sequence_array_ptr = std::make_unique<TranscriptionSequenceArray>();

  for (auto const& [offset, transcript_ptr] : getMap()) {

    auto sequence_validity = TranscriptionSequence::checkSequenceStatus(transcript_ptr);
    if (TranscriptionSequence::checkValidSequence(sequence_validity)) {

      if (not sequence_array_ptr->insertSequence(transcript_ptr)) {

        ExecEnv::log().error("Duplicate offset sequence, gene id: {}, offset: {}", transcript_ptr->getGene()->id(), offset);

      }

    }

  }

  return sequence_array_ptr;

}



