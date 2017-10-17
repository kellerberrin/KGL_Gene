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
// FeatureRecord members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::FeatureRecord::addSuperFeature(const kgl::FeatureIdent_t &super_feature_id,
                                         const std::shared_ptr<kgl::FeatureRecord> &super_feature_ptr) {

  super_features_.insert(std::make_pair(super_feature_id, super_feature_ptr));

}


void kgl::FeatureRecord::addSubFeature(const FeatureIdent_t& sub_feature_id,
                                       const std::shared_ptr<FeatureRecord>& sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}

