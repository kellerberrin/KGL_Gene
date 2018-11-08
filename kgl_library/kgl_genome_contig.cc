//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

  return true;

}

bool kgl::ContigFeatures::findFeatureId(const FeatureIdent_t& feature_id,
                                        std::vector<std::shared_ptr<kgl::Feature>>& feature_ptr_vec) const {

  auto iter_pair = id_feature_map_.equal_range(feature_id);

  feature_ptr_vec.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    feature_ptr_vec.emplace_back(iter->second);

  }

  return not feature_ptr_vec.empty();

}



void kgl::ContigFeatures::setupFeatureHierarchy() {

  // Remove all hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    Feature& feature = *feature_pair.second;
    // Remove feature hierarchy.
    feature.clearHierachy();
  }

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    Feature& feature = *feature_pair.second;
    // For each feature lookup a list of super_features
    std::vector<FeatureIdent_t> super_features;
    feature.getAttributes().getSuperFeatureIds(super_features);
    // TSS features are assigned to GENES
    std::vector<FeatureIdent_t> assigned_features;
    feature.getAttributes().getAssignedFeatureIds(assigned_features);
    // Merge TSS assigned features into super features.
    super_features.insert( super_features.end(), assigned_features.begin(), assigned_features.end() );
    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : super_features) {

      std::vector<std::shared_ptr<kgl::Feature>> super_feature_ptr_vec;
      if (not findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // If the feature is a TSS and is unassigned then continue.
        if (feature.isTSS() and super_feature_id == TSSFeature::TSS_UNASSIGNED) {

          continue;

        }

        // Otherwise flag an Error; could not find super feature.
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



bool kgl::ContigFeatures::findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const {

  gene_ptr_vec.clear();
  auto lb_result = gene_map_.lower_bound(offset);

  if (lb_result == gene_map_.end()) {

    return false;

  }

  auto result = gene_map_.equal_range(lb_result->first);

  for (auto it = result.first; it != result.second; ++it) {

    if (offset >= it->second->sequence().begin()) {

      gene_ptr_vec.emplace_back(it->second);

    }

  }

  return not gene_ptr_vec.empty();

}


// Convenience routine for Amino sequences.
std::shared_ptr<kgl::AminoSequence>
kgl::ContigFeatures::getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const {

  return coding_table_.getAminoSequence(sequence_ptr);

}


// Given a gene id and an mRNA id (sequence id) return the coding base sequence.
bool kgl::ContigFeatures::getCodingSequence(const FeatureIdent_t& gene_id,
                                            const FeatureIdent_t& sequence_id,
                                            std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const {

  std::vector<std::shared_ptr<Feature>> feature_ptr_vec;
  std::shared_ptr<const GeneFeature> gene_ptr;
  if (findFeatureId(gene_id, feature_ptr_vec)) {

    for (const auto& feature_ptr : feature_ptr_vec) {

      if (feature_ptr->id() == gene_id) {

        if (feature_ptr->isGene()) {

          gene_ptr = std::dynamic_pointer_cast<const GeneFeature>(feature_ptr);
          break;

        } else {

          ExecEnv::log().warn("Feature: {} is not a gene.", gene_id);
          return false;

        }

      }

    }

  }

  if (not gene_ptr) {

    ExecEnv::log().warn("Gene not found for feature id: {}.", gene_id);
    return false;

  }

  std::shared_ptr<const CodingSequenceArray> sequence_array_ptr = GeneFeature::getCodingSequences(gene_ptr);

  if (sequence_array_ptr->empty()) {

    ExecEnv::log().warn("No valid coding sequences found for Gene: {}.", gene_id);
    return false;

  }

  for (const auto& sequence: sequence_array_ptr->getMap()) {

    if (sequence.second->getCDSParent()->id() == sequence_id) {

      coding_sequence_ptr = sequence.second;
      return true;

    }

  }

  ExecEnv::log().warn("No valid coding sequences found for sequence id: {}.", sequence_id);
  return false;

}


// Given a CDS coding sequence, return the corresponding DNA base sequence (strand adjusted).
bool kgl::ContigFeatures::getDNA5SequenceCoding(const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                                                std::shared_ptr<DNA5SequenceCoding>& sequence_ptr) const {

  if (coding_sequence_ptr) {

    sequence_ptr = sequence_ptr_->DNA5SequenceContig::codingSequence(coding_sequence_ptr);
    return true;

  }

  ExecEnv::log().error("getDNA5SequenceCoding(), coding_sequence_ptr is null");
  return false;

}


kgl::TSSVector kgl::ContigFeatures::getTSSVector() const {

  TSSVector tss_vector;

  for (auto feature : offset_feature_map_) {

    if (not feature.second->isTSS()) continue;

    std::shared_ptr<const TSSFeature> tss_feature = std::dynamic_pointer_cast<TSSFeature>(feature.second);

    if (not tss_feature) {

      ExecEnv::log().error("Unexpected feature type for feature: {}", feature.second->id());

    } else {

      tss_vector.push_back(tss_feature);

    }

  }

  return tss_vector;

}
