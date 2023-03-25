//
// Created by kellerberrin on 3/12/18.
//

#include "kgl_genome_contig_aux.h"
#include "kgl_resource_db.h"
#include "kel_patterns.h"


namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AdjalleyTSSFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::AdjalleyTSSFeatures::checkAddFeature(std::shared_ptr<Feature>& feature_ptr) {

  if (not feature_ptr->isTSS()) {

    ExecEnv::log().error("AdjalleyTSSFeatures::checkAddFeature; Feature: {}, can only add TSS features to the AdjalleyTSSFeatures structure",
                         feature_ptr->id());
    return false;

  }

  addFeature(feature_ptr);

  return true;

}


void kgl::AdjalleyTSSFeatures::setupVerifyHierarchy(const StructuredFeatures& gene_super_features) {

  setupFeatureHierarchy(gene_super_features);
  verifyFeatureHierarchy();

}


std::vector<std::shared_ptr<const kgl::Feature>> kgl::AdjalleyTSSFeatures::getTSSVector() const {

  std::vector<std::shared_ptr<const kgl::Feature>> tss_vector;

  for (auto [offset, feature_ptr] : offsetFeatureMap()) {

    if (feature_ptr->isTSS()) {

      tss_vector.push_back(feature_ptr);

    } else {

      ExecEnv::log().error("Unexpected feature type: {} for TSS feature: {}", feature_ptr->featureType(), feature_ptr->id());

    }

  }

  return tss_vector;

}


void kgl::AdjalleyTSSFeatures::setupFeatureHierarchy(const StructuredFeatures& gene_super_features) {

  // Remove all hierarchies for all features.
  clearHierarchy();

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : idFeatureMap()) {

    Feature& feature = *feature_pair.second;

    // TSS features are assigned to GENES
    std::vector<FeatureIdent_t> assigned_features;
    feature.getAttributes().getAssignedFeatureIds(assigned_features);

    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : assigned_features) {

      std::vector<std::shared_ptr<const Feature>> super_feature_ptr_vec;
      if (not gene_super_features.findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // Otherwise flag an Error; could not find super feature.
        ExecEnv::log().warn("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      if (not super_feature_ptr_vec.empty()) {

        feature.setSuperFeature(super_feature_ptr_vec.front());
        std::const_pointer_cast<Feature>(super_feature_ptr_vec.front())->addSubFeature(feature.id(), feature_pair.second);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A collection of objects that contain auxillary genome information sources e.g. TSS, Histone modification, Promoter sites etc.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::AuxContigFeatures::setupVerifyHierarchy(const StructuredFeatures& gene_super_features) {

  adjalley_TSS_Features_.setupVerifyHierarchy(gene_super_features);

}


void kgl::AuxContigFeatures::checkAddFeature(std::shared_ptr<Feature>& feature_ptr) {

  if (not adjalley_TSS_Features_.checkAddFeature(feature_ptr)) {

    ExecEnv::log().error("AuxContigFeatures::checkAddFeature(), Problem adding TSS feature: {}", feature_ptr->id());

  }

}