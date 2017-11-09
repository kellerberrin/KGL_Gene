//
// Created by kellerberrin on 9/11/17.
//

#include <algorithm>
#include "kgl_attributes.h"

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Attributes members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns false if key not found.
bool kgl::Attributes::getAttributes(const std::string& key, std::vector<std::string>& values) const {

  auto iter_pair = attributes_.equal_range(key);

  values.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    values.emplace_back(iter->second);

  }

  return not values.empty();

}


// Always succeeds; keys are uppercase.
void kgl::Attributes::insertAttribute(const std::string& key, const std::string& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  std::string upper_case_key = key;
  std::transform(upper_case_key.begin(), upper_case_key.end(), upper_case_key.begin(), ::toupper);
  attributes_.insert(std::make_pair(upper_case_key, value));

}


void kgl::Attributes::getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const {

  all_key_value_pairs.clear();

  for (auto key_value_pair : attributes_) {

    all_key_value_pairs.emplace_back(key_value_pair);

  }

}

