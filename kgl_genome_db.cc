// MIT License
//
// Copyright (c) 2017
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// Created by kellerberrin on 7/10/17.
//

#include "kgl_genome_db.h"
#include "kgl_exec_env.h"

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
// FeatureRecord members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::FeatureRecord::addParent(const kgl::FeatureIdent_t& parent_id,
                                   const std::shared_ptr<kgl::FeatureRecord>& parent_ptr) {

  parents_.insert(std::make_pair(parent_id, parent_ptr));

}


void kgl::FeatureRecord::addSubFeature(const FeatureIdent_t& sub_feature_id,
                                       const std::shared_ptr<FeatureRecord>& sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigRecord members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigRecord::addFeature(std::shared_ptr<kgl::FeatureRecord>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

  return true;

}

bool kgl::ContigRecord::findFeatureId(kgl::FeatureIdent_t& feature_id,
                                      std::vector<std::shared_ptr<kgl::FeatureRecord>>& feature_ptr_vec) {

  auto iter_pair = id_feature_map_.equal_range(feature_id);

  feature_ptr_vec.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    feature_ptr_vec.emplace_back(iter->second);

  }

  return not feature_ptr_vec.empty();

}



void kgl::ContigRecord::setupFeatureHierarchy() {

  // Remove all hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord& feature = *feature_pair.second;
    // Remove feature hierarchy.
    feature.clearHierachy();
  }

  // Re-establish the hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord& feature = *feature_pair.second;
    // For each feature lookup a list of parents
    std::vector<FeatureIdent_t> parents;
    feature.getAttributes().getParents(parents);
    // Add parent pointers for the child and child pointers for the parents.
    for (auto parent_id : parents) {

      std::vector<std::shared_ptr<kgl::FeatureRecord>> parent_ptr_vec;
      if (not findFeatureId(parent_id, parent_ptr_vec)) {

        // Error; could not find parent.
        kgl::ExecEnv::log().error("Child Feature: {}; Parent Feature: {} does not exist", feature.id(), parent_id);

      }
      if (parent_ptr_vec.size() > 1) {

        // Warning, more than 1 parent.
        kgl::ExecEnv::log().warn("Parent id: {} returned : {} Parent Features", parent_id, parent_ptr_vec.size());

      }
      for (auto& parent_ptr : parent_ptr_vec) {

        feature.addParent(parent_id, parent_ptr);
        parent_ptr->addSubFeature(feature.id(), feature_pair.second);

      } // For all parents with same id.

    } // For all parent ids.

  } // For all features.

}


void kgl::ContigRecord::verifyFeatureHierarchy() {

  // No features larger than the contig.
  ContigSize_t contig_size = contigSize();
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord& feature = *feature_pair.second;
    // Error if feature overlaps the and of the contig.
    if (feature.sequence().end() >= contig_size) {

      kgl::ExecEnv::log().error("Feature: {}; (zero offset) end sequence: {} >= contig size: {}",
                                feature.id(), feature.sequence().end(), contig_size);

    }

  }

  // Check that children fit within parents
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord& feature = *feature_pair.second;

    for(auto child_feature_pair : feature.subFeatures()) {
      FeatureRecord &child_feature = *child_feature_pair.second;
      if (child_feature.sequence().begin() < feature.sequence().begin()
          or child_feature.sequence().end() > feature.sequence().end()) {

        kgl::ExecEnv::log().error("Child Feature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  child_feature.id(),
                                  child_feature.sequence().begin(),
                                  child_feature.sequence().end(),
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end());


      }

    }

    for(auto parent_feature_pair : feature.parentFeatures()) {
      FeatureRecord& parent_feature = *parent_feature_pair.second;
      if (feature.sequence().begin() < parent_feature.sequence().begin()
          or feature.sequence().end() > parent_feature.sequence().end()) {

        kgl::ExecEnv::log().error("Feature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end(),
                                  parent_feature.id(),
                                  parent_feature.sequence().begin(),
                                  parent_feature.sequence().end());

      } // If parent >= feature

    } // for all parents

    long duplicates = 0;
    // Check for duplicate parents.
    for(auto parent_feature_pair_a : feature.parentFeatures()) {

      for(auto parent_feature_pair_b : feature.parentFeatures()) {

        if (parent_feature_pair_a == parent_feature_pair_b) ++duplicates;

      } // for all parents

    } // for all parents

    duplicates -= feature.parentFeatures().size();

    if (duplicates > 0) {

      kgl::ExecEnv::log().info("Feature: {}; has {} duplicate parents", feature.id(), duplicates);

    }

    duplicates = 0;
    // Check for duplicate children.
    for(auto child_feature_pair_a : feature.subFeatures()) {

      for(auto child_feature_pair_b : feature.subFeatures()) {

        if (child_feature_pair_a == child_feature_pair_b) ++duplicates;

      } // for all children

    } // for all children

    duplicates -= feature.subFeatures().size();

    if (duplicates > 0) {

      kgl::ExecEnv::log().info("Feature: {}; has {} duplicate children", feature.id(), duplicates);

    }

  } // for all features.

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeSequences members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GenomeSequences::addContigSequence(const kgl::ContigId_t& contig_id, kgl::Sequence_t sequence) {

  using ContigPtr = std::shared_ptr<kgl::ContigRecord>;
  ContigPtr contig_ptr(std::make_shared<kgl::ContigRecord>(contig_id, std::move(sequence)));

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));

  return result.second;

}

bool kgl::GenomeSequences::getContigSequence(const kgl::ContigId_t& contig_id,
                                             std::shared_ptr<ContigRecord>& contig_ptr) const {

  auto result_iter = genome_sequence_map_.find(contig_id);

  if (result_iter != genome_sequence_map_.end()) {

    contig_ptr = result_iter->second;
    return true;

  }

  return false;

}

void kgl::GenomeSequences::createVerifyGenomeDatabase() {

  setupFeatureHierarchy();
  verifyFeatureHierarchy();

}

void kgl::GenomeSequences::setupFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setupFeatureHierarchy();

  }

}


void kgl::GenomeSequences::verifyFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyFeatureHierarchy();

  }

}

