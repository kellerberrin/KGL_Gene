//
// Created by kellerberrin on 10/10/17.
//

#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"
#include "kgl_alphabet_amino.h"

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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Feature::addSuperFeature(const kgl::FeatureIdent_t &super_feature_id,
                                         const std::shared_ptr<kgl::Feature> &super_feature_ptr) {

  super_features_.insert(std::make_pair(super_feature_id, super_feature_ptr));

}


void kgl::Feature::addSubFeature(const FeatureIdent_t& sub_feature_id,
                                       const std::shared_ptr<Feature>& sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}


void kgl::Feature::recusivelyPrintsubfeatures(long feature_level) const {

  for (const auto& feature : sub_features_) {

    ExecEnv::log().info("Level: {}, Feature: {}, Type: {}, begin: {}, end: {} strand: {}",
                        feature_level,
                        feature.second->id(),
                        feature.second->featureType(),
                        feature.second->sequence().begin(),
                        feature.second->sequence().end(),
                        static_cast<char>(feature.second->sequence().strand()));

    feature.second->recusivelyPrintsubfeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}



bool kgl::Feature::verifyCDSPhase(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const {

  bool result = true;
  // Check for mod3
  for(const auto& sorted_cds : coding_seq_ptr->getMap()) {

    result = result and verifyMod3(sorted_cds.second->getSortedCDS());
    result = result and verifyStrand(sorted_cds.second->getSortedCDS());
    result = result and verifyPhase(sorted_cds.second->getSortedCDS());

  }

  return result;

}


bool kgl::Feature::verifyMod3(const SortedCDS& sorted_cds) const {

  bool result = true;
// Check the combined sequence length is mod 3 = 0

  ContigSize_t coding_sequence_length = 0;
  for (auto cds : sorted_cds) {

    coding_sequence_length += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  if ((coding_sequence_length % AminoAcidTypes::CODON_SIZE) != 0) {

    ExecEnv::log().warn("Gene: {} offset: {} CDS coding sequence length mod 3 not zero : {}",
                        id(),
                        sequence().begin(),
                        (coding_sequence_length % 3));

    result = false;

  }

  return result;

}

bool kgl::Feature::verifyStrand(const SortedCDS& sorted_cds) const {

  bool result = true;

// Check the strand is consistent and not unknown.
  for (auto cds : sorted_cds) {

    if (cds.second->sequence().strand() != sequence().strand()) {

      ExecEnv::log().error("CDS: {} offset: {} strand: {}, parent sequence strand: {} mis-match",
                           cds.second->id(),
                           cds.second->sequence().begin(),
                           static_cast<char>(cds.second->sequence().strand()),
                           static_cast<char>(sequence().strand()));
      result = false;

    }

  }

  return result;

}


bool kgl::Feature::verifyPhase(const SortedCDS& sorted_cds) const {

  bool result = true;

  switch(sequence().strand()) {

    case StrandSense::FORWARD: {

      ContigOffset_t sequence_length = 0;
      bool warn_adjust = false;
      static bool warn_forward_once = false;
      for (auto it = sorted_cds.begin(); it != sorted_cds.end(); ++it) {

        CDSPhaseType_t phase = (AminoAcidTypes::CODON_SIZE - (sequence_length % AminoAcidTypes::CODON_SIZE))
                               % AminoAcidTypes::CODON_SIZE;

        if (it->second->phase() != phase) {

          std::const_pointer_cast<CDSFeature>(it->second)->phase(phase);
          warn_adjust = true;

        }

        sequence_length += it->second->sequence().end() - it->second->sequence().begin();

      }

      if (warn_adjust and not warn_forward_once) {
        warn_forward_once = true;
        ExecEnv::log().warn("Inconsistent Gff3 phase fields found in forward ('+') strand CDS features - adjusted");
      }
    }
    break;

    case StrandSense::REVERSE: {

      ContigOffset_t sequence_length = 0;
      bool warn_adjust = false;
      static bool warn_reverse_once = false;
      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        CDSPhaseType_t phase = (AminoAcidTypes::CODON_SIZE - (sequence_length % AminoAcidTypes::CODON_SIZE))
                               % AminoAcidTypes::CODON_SIZE;

        if (rit->second->phase() != phase) {

          std::const_pointer_cast<CDSFeature>(rit->second)->phase(phase);
          warn_adjust = true;

        }

        sequence_length += rit->second->sequence().end() - rit->second->sequence().begin();

      }

      if (warn_adjust and not warn_reverse_once) {

        warn_reverse_once = true;
        ExecEnv::log().warn("Inconsistent Gff3 phase fields found in reverse ('-') strand CDS features - adjusted");

      }

    }
    break;

    case StrandSense::UNKNOWN:
    ExecEnv::log().error("verifyPhase() CDS parent: {} sequence strand is UNKNOWN '.'", id());
    result = false;
    break;

  }

  return result;

}


std::shared_ptr<kgl::Feature> kgl::Feature::getGene() const {

  // recursivly search upward

  for (auto feature : super_features_) {

    if (feature.second->isGene()) {

      return feature.second;

    } else {

      return feature.second->getGene();

    }
  }

  std::shared_ptr<kgl::Feature> null;
  return null;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<const kgl::CodingSequenceArray>
kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene) {

  std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>());
  getCodingSequences(gene, gene, sequence_array_ptr);
  return sequence_array_ptr;

}

// This routine is recursive. Assumes all the CDS are on the same sub-feature level.
bool kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene_ptr,
                                          std::shared_ptr<const Feature> cds_parent_ptr,
                                          std::shared_ptr<CodingSequenceArray>& sequence_array_ptr) {

  bool result = true;
  SortedCDS parent_cds;

  for (auto sub_feature : cds_parent_ptr->subFeatures()) {

    if (sub_feature.second->isCDS()) {


      auto insert = parent_cds.insert(std::make_pair(sub_feature.second->sequence().begin(),
                                                     std::static_pointer_cast<const CDSFeature>(sub_feature.second)));

      if (not insert.second) {

        ExecEnv::log().warn("Duplicate CDS: {} at contig offset: {}",
                            sub_feature.second->id(),
                            sub_feature.second->sequence().begin());
        result = false;
      }

    } else { // recursively call this function

      result = result and getCodingSequences(gene_ptr, sub_feature.second, sequence_array_ptr);

    }

  }

  if (not parent_cds.empty()) {

    std::shared_ptr<const CodingSequence> coding_sequence(std::make_shared<const CodingSequence>(gene_ptr,
                                                                                                 cds_parent_ptr,
                                                                                                 parent_cds));

    result = result and sequence_array_ptr->insertCodingSequence(coding_sequence);

  }

  return result;

}


std::shared_ptr<const kgl::CodingSequenceArray>
kgl::GeneFeature::getOffsetSequences(ContigOffset_t offset,
                                     std::shared_ptr<const CodingSequenceArray> sequence_array_ptr) {

  std::shared_ptr<CodingSequenceArray> filtered_map_ptr(std::make_shared<CodingSequenceArray>());
  for (const auto& sequence : sequence_array_ptr->getMap()) {

    if (sequence.second->isWithinCoding(offset)) {

      filtered_map_ptr->insertCodingSequence(sequence.second);

    }

  }

  return filtered_map_ptr;

}

