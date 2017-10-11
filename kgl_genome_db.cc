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

#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


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

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord& feature = *feature_pair.second;
    // For each feature lookup a list of super_features
    std::vector<FeatureIdent_t> super_features;
    feature.getAttributes().getSuperFeatureIds(super_features);
    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : super_features) {

      std::vector<std::shared_ptr<kgl::FeatureRecord>> super_feature_ptr_vec;
      if (not findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // Error; could not find super feature.
        kgl::ExecEnv::log().error("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        kgl::ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      for (auto& super_feature_ptr : super_feature_ptr_vec) {

        feature.addSuperFeature(super_feature_id, super_feature_ptr);
        super_feature_ptr->addSubFeature(feature.id(), feature_pair.second);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

}


void kgl::ContigRecord::verifyFeatureHierarchy() {

  verifyContigOverlap();
  verifySubFeatureSuperFeatureDimensions();
  removeSubFeatureDuplicates();
  removeSuperFeatureDuplicates();
  verifySubFeatureDuplicates();
  verifySuperFeatureDuplicates();

}


void kgl::ContigRecord::verifyContigOverlap() {

  // If feature dimensions are [1, size] instead of [0, size-1] then assume that conversion from the
  // Gff convention of [1, size] has not been performed correctly during feature read from disk.
  // Adjust to [0, size-1] here.
  // Note that this suggests a problem with the (3rd party) Gff read functionality and should be addressed there.

  long adjust_from_1_size = 0;  // Report heuristic adjustment.

  // No features larger than the contig.
  ContigSize_t contig_size = contigSize();

  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;
    // Error if feature overlaps the and of the contig.
    // If [1,contig_size] then adjust to [0, contis_size-1]
    if (feature.sequence().end() >= contig_size) {

      if (feature.sequence().end() == contig_size and feature.sequence().begin() <= 1) { // adjust to [0, size-1]

        ++adjust_from_1_size;
        FeatureSequence adj_sequence = feature.sequence();
        adj_sequence.end(adj_sequence.end() - 1);
        adj_sequence.begin(0);
        feature.sequence(adj_sequence);

      } else {

        kgl::ExecEnv::log().error("Feature: {};  sequence: [{}:{}] >= contig: {} size: {}",
                                  feature.id(), feature.sequence().begin(), feature.sequence().end(),
                                  contigId(), contig_size);

      } // if end == contig_size

    } // contig overlap

  } // for all features.

  if (adjust_from_1_size > 0) {

    kgl::ExecEnv::log().warn("Contig: {}; size: {} had: {} 1-offset features [1, {}], adjusted to zero-offset [0, {}]",
                             contigId(), contig_size, adjust_from_1_size, contig_size, (contig_size - 1));
  }
}


void kgl::ContigRecord::verifySubFeatureSuperFeatureDimensions() {


  // Check that sub-features fit within feature.
  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;

    for (auto sub_feature_pair : feature.subFeatures()) {
      FeatureRecord &sub_feature = *sub_feature_pair.second;
      if (sub_feature.sequence().begin() < feature.sequence().begin()
          or sub_feature.sequence().end() > feature.sequence().end()) {

        kgl::ExecEnv::log().error("SubFeature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  sub_feature.id(),
                                  sub_feature.sequence().begin(),
                                  sub_feature.sequence().end(),
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end());


      } // if sub-feature overlaps feature

    } // for all sub-features.

    // Check that features fit within super-feature.
    for (auto super_feature_pair : feature.superFeatures()) {
      FeatureRecord &super_feature = *super_feature_pair.second;
      if (feature.sequence().begin() < super_feature.sequence().begin()
          or feature.sequence().end() > super_feature.sequence().end()) {

        kgl::ExecEnv::log().error("Feature: {}; [{}:{}] overlaps SuperFeature {}; [{}:{}]",
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end(),
                                  super_feature.id(),
                                  super_feature.sequence().begin(),
                                  super_feature.sequence().end());

      } // If feature overlaps super-feature.

    } // for all super features.

  } // for all features.

}


void kgl::ContigRecord::removeSubFeatureDuplicates() {

  long duplicates_removed = 0;

  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;

    duplicates_removed += deleteIterable(feature.subFeatures());

  }

  if (duplicates_removed > 0) {

    kgl::ExecEnv::log().info("{} duplicate sub-features removed from contig: {}", duplicates_removed, contigId());

  }

}

void kgl::ContigRecord::removeSuperFeatureDuplicates() {

  long duplicates_removed = 0;

  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;

    duplicates_removed += deleteIterable(feature.superFeatures());

  }

  if (duplicates_removed > 0) {

    kgl::ExecEnv::log().info("{} duplicate super-features removed from contig: {}", duplicates_removed, contigId());

  }

}


void kgl::ContigRecord::verifySubFeatureDuplicates() {

  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.subFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate sub-features", feature.id(), duplicates);

    }

  }

}


void kgl::ContigRecord::verifySuperFeatureDuplicates() {

  for (auto feature_pair : id_feature_map_) {
    FeatureRecord &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.superFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate super-features", feature.id(), duplicates);

    }

  }

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


void kgl::GenomeSequences::registerContigData(std::shared_ptr<kgl::ContigDataBlock>& contig_data_ptr) {
// Create data blocks for each contig in the genome database
  for (const auto &contig_pair : getGenomeSequenceMap()) {

    if (not contig_data_ptr->insertContig(contig_pair.first, contig_pair.second->sequence().length())) {

      kgl::ExecEnv::log().error("ContigDataBlock; attempted to add duplicate contig; {}", contig_pair.first);

    }

  }

}