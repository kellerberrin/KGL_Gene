//
// Created by kellerberrin on 7/11/17.
//

#include "kel_exec_env.h"
#include "kgl_genome_feature.h"

namespace kgl = kellerberrin::genome;

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
// CodingSequence - Members
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::StrandSense kgl::CodingSequence::strand() const {

  return getGene()->sequence().strand();

}




// Offset of the start of the sequence - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::CodingSequence::start() const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().begin();
  }

  return coding_feature_map_.begin()->second->sequence().begin();

}

// Offset of the end of the sequence (last nucleotide + 1) - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::CodingSequence::end() const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().end();
  }

  return coding_feature_map_.rbegin()->second->sequence().end();

}


// offset of the nucleotide not in the coding sequence (strand adjusted).
kgl::ContigOffset_t kgl::CodingSequence::prime_5() const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().strand() == StrandSense::REVERSE ? getGene()->sequence().end() : getGene()->sequence().begin() - 1;
    return 0;
  }

  switch(strand()) {

    case StrandSense::FORWARD:
      return coding_feature_map_.begin()->second->sequence().begin() - 1;

    case StrandSense::REVERSE:
      return coding_feature_map_.rbegin()->second->sequence().end();

  }

  return coding_feature_map_.begin()->second->sequence().begin(); // never reached; to keep the compiler happy

}



void kgl::CodingSequence::prime_5_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    begin_offset = 0;
    size = 0;
    return;

  }

  switch(strand()) {

    case StrandSense::FORWARD:
      begin_offset = coding_feature_map_.begin()->second->sequence().begin() - requested_size;
      size = requested_size;
      break;

    case StrandSense::REVERSE:
      begin_offset = coding_feature_map_.rbegin()->second->sequence().end();
      size = requested_size;
      break;

  }


}



// offset of the nucleotide not in the coding seuence (strand adjusted).
kgl::ContigOffset_t kgl::CodingSequence::prime_3() const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("CodingSequence::prime_3; coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().strand() == StrandSense::REVERSE ? getGene()->sequence().begin() - 1 : getGene()->sequence().end();
  }


  switch(strand()) {

    case StrandSense::FORWARD:
      return coding_feature_map_.rbegin()->second->sequence().end();

    case StrandSense::REVERSE:
      return coding_feature_map_.begin()->second->sequence().begin() - 1;

  }

  return coding_feature_map_.begin()->second->sequence().begin(); // never reached; to keep the compiler happy

}


void kgl::CodingSequence::prime_3_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const {

  // Safety first.
  if (coding_feature_map_.empty()) {

    ExecEnv::log().error("CodingSequence::prime_3_region; coding sequence for gene id: {} is empty", getGene()->id());
    begin_offset = 0;
    size = 0;
    return;

  }


  switch(strand()) {

    case StrandSense::FORWARD:
      begin_offset = coding_feature_map_.rbegin()->second->sequence().end();
      size = requested_size;
      break;

    case StrandSense::REVERSE:
      begin_offset = coding_feature_map_.begin()->second->sequence().begin() - requested_size;
      size = requested_size;
      break;

  }


}


std::shared_ptr<const kgl::ContigReference> kgl::CodingSequence::contig() const {

  return getGene()->contig();

}



kgl::ContigSize_t kgl::CodingSequence::codingNucleotides() const {

  ContigSize_t coding_size = 0;
  // Loop through and test membership of each cds. Reminder; testing for [begin, end)
  for (const auto& cds : coding_feature_map_) {

    coding_size += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  return coding_size;

}



void kgl::CodingSequenceArray::printSequence(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) {

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
kgl::TranscriptionSequenceType kgl::CodingSequence::codingType() const {

  if (coding_feature_map_.empty()) {

    return TranscriptionSequenceType::EMPTY;

  }

  auto const& [offset, feature_ptr]  = *coding_feature_map_.begin();

  if (feature_ptr->superType() == Feature::CDS_TYPE_) {

    return TranscriptionSequenceType::PROTEIN;

  } else {

    return TranscriptionSequenceType::NCRNA;

  }

}

bool kgl::CodingSequenceArray::insertCodingSequence(std::shared_ptr<const CodingSequence> coding_sequence_ptr) {

#ifdef CODING_SEQUENCE_ISMAP // Using a map

  auto insert = coding_sequence_map_.insert(std::make_pair(coding_sequence_ptr->getCDSParent()->id(),
                                                           coding_sequence_ptr));

  if (not insert.second) {

    ExecEnv::log().warn("Duplicate CDS parent: {} at contig offset: {}",
                        coding_sequence_ptr->getCDSParent()->id(),
                        coding_sequence_ptr->getCDSParent()->sequence().begin());
  }

  return insert.second;

#else // Using a multimap. Same parent features permitted 

  coding_sequence_map_.insert(std::make_pair(coding_sequence_ptr->getParent()->id(), coding_sequence_ptr));

  return true;

#endif

}

// Assumes that all coding sequences are the same type.
kgl::TranscriptionSequenceType kgl::CodingSequenceArray::codingType() const {

  if (coding_sequence_map_.empty()) {

    return TranscriptionSequenceType::EMPTY;

  }

  auto const& first_sequence = getFirst();

  return first_sequence->codingType();

}


std::shared_ptr<const kgl::CodingSequence> kgl::CodingSequenceArray::getFirst() const {

  if (empty()) {

    ExecEnv::log().critical("CodingSequenceArray::getFirst() tried to obtain a CodingSequence from empty container");

  }

  return getMap().begin()->second;

}

