//
// Created by kellerberrin on 10/10/17.
//

#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

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

// Recursively descends the sub-feature tree and returns a vector of sorted
// CDS trees. Each sorted CDS tree is an alternative coding sequence for the gene.
// The logic in this function must take into account that the CDS features
// can be direct sub-features of the Gene or are sub-features of multiple mRna
// sub-features. This corresponds to the Gff3 standard.
// All this makes the logic a bit tricky, so read the code carefully before modifying.
bool kgl::Feature::getSortedCDS(SortedCDSVector& sorted_cds_vec) const {

  bool result = true;

  auto sub_feature_iter = sub_features_.begin();
  while (sub_feature_iter != sub_features_.end()) {

    if (sub_feature_iter->second->featureType() == CDSFeature::CDS_TYPE) {

      SortedCDS sorted_cds;

      while(sub_feature_iter != sub_features_.end()) {

        if (sub_feature_iter->second->featureType() == CDSFeature::CDS_TYPE) {

          auto insert_result = sorted_cds.insert(std::make_pair(sub_feature_iter->second->sequence().begin(),
                                                  std::static_pointer_cast<CDSFeature>(sub_feature_iter->second)));

          if (not insert_result.second) {

            ExecEnv::log().warn("Duplicate CDS: {} at contig offset: {}",
                                sub_feature_iter->second->id(),
                                sub_feature_iter->second->sequence().begin());
            result = false;
          }


        } // If sub-feature is a CDS

        ++sub_feature_iter; // next sub-feature.

      } // While looking for other CDS for this nRNA or Gene.

      if (not sorted_cds.empty()) { // Unnecessary, but safe.

        sorted_cds_vec.emplace_back(sorted_cds);

      }

    } else { // Not a CDS, so look deeper

      sub_feature_iter->second->getSortedCDS(sorted_cds_vec); // Recursive call for the sub-feature.
      ++sub_feature_iter; // while next sub-feature.

    }

  }

  return result;

}


void kgl::Feature::recusivelyPrintsubfeatures(long feature_level) const {

  for (const auto& feature : sub_features_) {

    ExecEnv::log().info("Level: {}, Feature: {}, Type: {}, begin: {}, end: {} strand: {}",
                        feature_level,
                        feature.second->id(),
                        feature.second->featureType(),
                        feature.second->sequence().begin(),
                        feature.second->sequence().end(),
                        static_cast<char>(feature.second->sequence().sense()));

    feature.second->recusivelyPrintsubfeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}


bool kgl::Feature::verifyCDSPhase(const SortedCDSVector& parent_sorted_vec) {

  bool result = true;
  // Check for mod3
  for(auto& sorted_cds : parent_sorted_vec) {

    result = result and verifyMod3(sorted_cds);
    result = result and verifyPhase(sorted_cds);

  }

  return result;

}


bool kgl::Feature::verifyMod3(const SortedCDS& sorted_cds) {

  bool result = true;
// Check the combined sequence length is mod 3 = 0

  ContigSize_t coding_sequence_length = 0;
  for (auto cds : sorted_cds) {

    coding_sequence_length += (cds.second->sequence().end() - cds.second->sequence().begin());

  }

  if ((coding_sequence_length % 3) != 0) {

    ExecEnv::log().warn("Gene: {} offset: {} CDS coding sequence length mod 3 not zero : {}",
                        id(),
                        sequence().begin(),
                        (coding_sequence_length % 3));

    result = false;

  }

  return result;

}

bool kgl::Feature::verifyPhase(const SortedCDS& sorted_cds) {

  bool result = true;

// Check the strand is consistent and not unknown.
  for (auto cds : sorted_cds) {

    if (cds.second->sequence().sense() != sequence().sense()) {

      ExecEnv::log().error("CDS: {} offset: {} strand: {}, parent sequence strand: {} mis-match",
                           cds.second->id(),
                           cds.second->sequence().begin(),
                           static_cast<char>(cds.second->sequence().sense()),
                           static_cast<char>(sequence().sense()));
      result = false;

    }

  }

  switch(sequence().sense()) {

    case StrandSense::FORWARD: {

        CDSPhaseType_t phase = 0;
        for (auto it = sorted_cds.begin(); it != sorted_cds.end(); ++it) {

          if (it->second->phase() != phase) {

            ExecEnv::log().error("CDS: {} , offset: {}, mis-match; calculated phase: {}, CDS phase: {}",
                                 it->second->id(),
                                 it->second->sequence().begin(),
                                 phase,
                                 it->second->phase());
            result = false;

          }
          phase = (3 - ((it->second->sequence().end() - it->second->sequence().begin()) % 3)) % 3;
        }
      }
      break;

    case StrandSense::REVERSE:
      break;

    case StrandSense::UNKNOWN:
      ExecEnv::log().error("verifyPhase() CDS parent sequence strand is UNKNOWN");
      result = false;
      break;

  }

  return result;

}
