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

// Always succeeds; keys are uppercase.
void kgl::Attributes::insertAttribute(std::string&& key, std::string&& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  attributes_.emplace(Utility::toupper(Utility::trimEndWhiteSpace(key)), value);

}

std::string kgl::Attributes::getHGNC() const {

  std::string hgnc_id;

  for (auto const& [key, attrib] : attributes_) {

    if (key == DBXREF_ and attrib.find(HGNC_) == 0) {

      hgnc_id = attrib.substr(std::string(HGNC_).length(), std::string::npos);
      hgnc_id = Utility::trimEndWhiteSpace(hgnc_id);

    }

  }

  return hgnc_id;

}