//
// Created by kellerberrin on 9/11/17.
//

#include "kgl_genome_attributes.h"
#include "kel_utility.h"
#include "kel_exec_env.h"

#include <algorithm>

namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Attributes members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns false if key not found.
bool kgl::Attributes::getAttributes(const std::string &key, std::vector<std::string> &values) const {

  auto iter_pair = attributes_.equal_range(key);

  values.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    values.push_back(iter->second);

  }

  return not values.empty();

}


// Always succeeds; keys are uppercase.
void kgl::Attributes::insertAttribute(const std::string& key, const std::string& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  attributes_.emplace(Utility::toupper(Utility::trimEndWhiteSpace(key)), value);

}


void kgl::Attributes::getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const {

  all_key_value_pairs.clear();

  for (auto key_value_pair : attributes_) {

    all_key_value_pairs.emplace_back(key_value_pair);

  }

}

