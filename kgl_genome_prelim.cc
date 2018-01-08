//
// Created by kellerberrin on 7/11/17.
//

#include "kgl_exec_env.h"
#include "kgl_genome_feature.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingSequence - Members
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::StrandSense kgl::CodingSequence::strand() const {

  return getGene()->sequence().strand();

}

// Offset of the start of the sequence - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::CodingSequence::start() const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().begin();
  }

  return sorted_cds_.begin()->second->sequence().begin();

}

// Offset of the end of the sequence (last nucleotide + 1) - not strand adjusted. Uses half interval [start, end).
kgl::ContigOffset_t kgl::CodingSequence::end() const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().end();
  }

  return sorted_cds_.rbegin()->second->sequence().end();

}


// offset of the nucleotide not in the coding sequence (strand adjusted).
kgl::ContigOffset_t kgl::CodingSequence::prime_5() const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().strand() == StrandSense::REVERSE ? getGene()->sequence().end() : getGene()->sequence().begin() - 1;
    return 0;
  }

  switch(strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().error("prime_5(), unknown strand sense for gene id: {}", getGene()->id());
    case StrandSense::FORWARD:
      return sorted_cds_.begin()->second->sequence().begin() - 1;

    case StrandSense::REVERSE:
      return sorted_cds_.rbegin()->second->sequence().end();

  }

  return sorted_cds_.begin()->second->sequence().begin(); // never reached; to keep the compiler happy

}



void kgl::CodingSequence::prime_5_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    begin_offset = 0;
    size = 0;
    return;

  }

  switch(strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().error("prime_5_region(), unknown strand sense for gene id: {}", getGene()->id());
    case StrandSense::FORWARD:
      begin_offset = sorted_cds_.begin()->second->sequence().begin() - requested_size;
      size = requested_size;
      break;

    case StrandSense::REVERSE:
      begin_offset = sorted_cds_.rbegin()->second->sequence().end();
      size = requested_size;
      break;

  }


}



// offset of the nucleotide not in the coding seuence (strand adjusted).
kgl::ContigOffset_t kgl::CodingSequence::prime_3() const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    return getGene()->sequence().strand() == StrandSense::REVERSE ? getGene()->sequence().begin() - 1 : getGene()->sequence().end();
  }

  switch(strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().error("prime_3(), unknown strand sense for gene id: {}", getGene()->id());
    case StrandSense::FORWARD:
      return sorted_cds_.rbegin()->second->sequence().end();

    case StrandSense::REVERSE:
      return sorted_cds_.begin()->second->sequence().begin() - 1;

  }

  return sorted_cds_.begin()->second->sequence().begin(); // never reached; to keep the compiler happy

}


void kgl::CodingSequence::prime_3_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const {

  // Safety first.
  if (sorted_cds_.empty()) {

    ExecEnv::log().error("prime_5(), coding sequence for gene id: {} is empty", getGene()->id());
    begin_offset = 0;
    size = 0;
    return;

  }

  switch(strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().error("prime_3_region(), unknown strand sense for gene id: {}", getGene()->id());
    case StrandSense::FORWARD:
      begin_offset = sorted_cds_.rbegin()->second->sequence().end();
      size = requested_size;
      break;

    case StrandSense::REVERSE:
      begin_offset = sorted_cds_.begin()->second->sequence().begin() - requested_size;
      size = requested_size;
      break;

  }


}


// return true if the contig_offset lies with a CDS.
bool kgl::CodingSequence::isWithinCoding(ContigOffset_t contig_offset) const {

  // Safety first.
  if (sorted_cds_.empty()) return false;

  // Less than the begin offset.
  if (contig_offset < sorted_cds_.begin()->second->sequence().begin()) return false;

  // More than the end offset - remember the end offset is 1 past the last nucleotide. [begin, end)
  if (contig_offset >= sorted_cds_.rbegin()->second->sequence().end()) return false;

  // Loop through and test membership of each cds. Reminder; testing for [begin, end)
  for (const auto& cds : sorted_cds_) {

    if (cds.second->sequence().begin() <= contig_offset and cds.second->sequence().end() > contig_offset) {

      return true;

    }

  }

  return false;

}


kgl::ContigSize_t kgl::CodingSequence::codingNucleotides() const {

  ContigSize_t coding_size = 0;
  // Loop through and test membership of each cds. Reminder; testing for [begin, end)
  for (const auto& cds : sorted_cds_) {

    coding_size += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  return coding_size;

}



void kgl::CodingSequenceArray::printCodingSequence(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) {

  long vector_count = 0;
  for (const auto& sequence : coding_seq_ptr->getMap()) {

    ExecEnv::log().info("Gene: {}, begin: {}, end: {} strand: {}",
                        sequence.second->getGene()->id(),
                        sequence.second->getGene()->sequence().begin(),
                        sequence.second->getGene()->sequence().end(),
                        static_cast<char>(sequence.second->getGene()->sequence().strand()));

    ExecEnv::log().info("Parent Feature: {}, begin: {}, end: {} strand: {}",
                        sequence.second->getCDSParent()->id(),
                        sequence.second->getCDSParent()->sequence().begin(),
                        sequence.second->getCDSParent()->sequence().end(),
                        static_cast<char>(sequence.second->getCDSParent()->sequence().strand()));

    ++vector_count;

    ExecEnv::log().info("++++++++++++++ CDS Vector : {} ********************", vector_count);

    for (const auto& cds : sequence.second->getSortedCDS()) {

      ExecEnv::log().info("CDS: {}, Type: {}, begin: {}, end: {} strand: {}",
                          cds.second->id(),
                          cds.second->featureType(),
                          cds.second->sequence().begin(),
                          cds.second->sequence().end(),
                          static_cast<char>(cds.second->sequence().strand()));

    }

  }

}


bool kgl::CodingSequenceArray::insertCodingSequence(std::shared_ptr<const CodingSequence> coding_sequence_ptr) {

  auto insert = coding_sequence_map_.insert(std::make_pair(coding_sequence_ptr->getCDSParent()->id(),
                                                           coding_sequence_ptr));
  if (not insert.second) {

    ExecEnv::log().warn("Duplicate CDS parent: {} at contig offset: {}",
                        coding_sequence_ptr->getCDSParent()->id(),
                        coding_sequence_ptr->getCDSParent()->sequence().begin());
  }

  return insert.second;

}


std::shared_ptr<const kgl::CodingSequence> kgl::CodingSequenceArray::getFirst() const {

  if (empty()) {

    ExecEnv::log().critical("CodingSequenceArray::getFirst() tried to obtain a CodingSequence from empty container");

  }

  return getMap().begin()->second;

}
