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


bool kgl::FeatureRecord::addParent(const kgl::FeatureIdent_t& parent_id,
                                   const std::shared_ptr<kgl::FeatureRecord>& parent_ptr) {

  auto result = sub_features_.insert(std::make_pair(parent_id, parent_ptr));

  return result.second;

}


bool kgl::FeatureRecord::addSubFeature(const FeatureIdent_t& sub_feature_id,
                                       const std::shared_ptr<FeatureRecord>& sub_feature_ptr) {

  auto result = sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

  return result.second;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigRecord members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigRecord::addFeature(std::shared_ptr<kgl::FeatureRecord>& feature_ptr) {

  auto result = id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  if (not result.second) {

    return false;

  }

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

  return true;

}


bool kgl::ContigRecord::findFeatureId(kgl::FeatureIdent_t& feature_id,
                                      std::shared_ptr<kgl::FeatureRecord>& feature_ptr) {

  auto result = id_feature_map_.find(feature_id);

  if (result != id_feature_map_.end()) {

    feature_ptr = result->second;
    return true;

  }

  return false;  // Not found.

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

      std::shared_ptr<kgl::FeatureRecord> parent_ptr;
      if (findFeatureId(parent_id, parent_ptr)) {

        if (not feature.addParent(parent_id, parent_ptr)) {

          // Error; parent already exists

        }

        if(not parent_ptr->addSubFeature(feature.id(), feature_pair.second)) {

          // Error; child already exists

        }

      } else {

        // Error; could not find parent.

      } // Find parent.

    } // For all parents.

  } // For all features.

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

void kgl::GenomeSequences::setupFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setupFeatureHierarchy();

  }

}
