//
// Created by kellerberrin on 1/12/18.
//


#include "kgl_genome_contig_feature.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// AuxGenomeFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool kgl::StructuredFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

  return true;

}

bool kgl::StructuredFeatures::findFeatureId(const FeatureIdent_t& feature_id,
                                           std::vector<std::shared_ptr<kgl::Feature>>& feature_ptr_vec) const {

  auto iter_pair = id_feature_map_.equal_range(feature_id);

  feature_ptr_vec.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    feature_ptr_vec.emplace_back(iter->second);

  }

  return not feature_ptr_vec.empty();

}



kgl::TSSVector kgl::TSSFeatures::getTSSVector() const {

  TSSVector tss_vector;

  for (auto feature : offetFeatureMap()) {

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

