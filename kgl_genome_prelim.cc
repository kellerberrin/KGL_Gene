//
// Created by kellerberrin on 7/11/17.
//

#include "kgl_exec_env.h"
#include "kgl_genome_feature.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureAttributes members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns false if key not found.
bool kgl::FeatureAttributes::getAttributes(const std::string& key, std::vector<std::string>& values) const {

  auto iter_pair = attributes_.equal_range(key);

  values.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    values.emplace_back(iter->second);

  }

  return not values.empty();

}


// Always succeeds; keys are uppercase.
void kgl::FeatureAttributes::insertAttribute(const std::string& key, const std::string& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  std::string upper_case_key = key;
  std::transform(upper_case_key.begin(), upper_case_key.end(), upper_case_key.begin(), ::toupper);
  attributes_.insert(std::make_pair(upper_case_key, value));

}


void kgl::FeatureAttributes::getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const {

  all_key_value_pairs.clear();

  for (auto key_value_pair : attributes_) {

    all_key_value_pairs.emplace_back(key_value_pair);

  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingSequence - Members
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// return true if the contig_offset lies with a CDS.
bool kgl::CodingSequence::isWithinCoding(ContigOffset_t contig_offset) const {

  // Safety first.
  if (sorted_cds_.empty()) return false;

  // Less than the begin offset.
  if (contig_offset < sorted_cds_.begin()->second->sequence().begin()) return false;

  // More than the end offset.
  if (contig_offset > sorted_cds_.rbegin()->second->sequence().end()) return false;

  // Loop through and test membership of each cds. Reminder; testing for [begin, end)
  for (const auto& cds : sorted_cds_) {

    if (cds.second->sequence().begin() <= contig_offset and cds.second->sequence().end() > contig_offset) {

      return true;

    }

  }

  return false;

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


void kgl::CodingSequenceArray::mergeArrays(std::shared_ptr<const CodingSequenceArray> merge_array) {

  for (auto sequence: merge_array->getMap()) {

    insertCodingSequence(sequence.second);

  }

}
