//
// Created by kellerberrin on 9/11/17.
//

#include "kgl_genome_attributes.h"
#include "kel_utility.h"
#include "kel_exec_env.h"

#include <ranges>
#include <iterator>
#include <string_view>


namespace kgl = kellerberrin::genome;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Attributes members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns false if key not found.
bool kgl::Attributes::getAttributes(const std::string &key, std::vector<std::string> &values) const {

  values.clear();
  auto const [first, last] = attributes_.equal_range(key);
  std::ranges::copy(std::ranges::subrange(first, last) | std::views::values, std::back_inserter(values));

  return not values.empty();

}


// Returns false if key not found.
std::vector<std::string> kgl::Attributes::getAttributes(const std::string &key) const {

  std::vector<std::string> values;
  getAttributes(key, values);
  return values;

}


// Always succeeds; keys are uppercased on insert but looked up verbatim.
void kgl::Attributes::insertAttribute(const std::string& key, const std::string& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  attributes_.emplace(Utility::toupper(Utility::trimEndWhiteSpace(key)), value);

}

// Always succeeds; keys are uppercased on insert but looked up verbatim.
void kgl::Attributes::insertAttribute(std::string&& key, std::string&& value) {

  // Convert the key to upper case to avoid the vagaries of non-standard case in keys.
  // The key must still be transformed, so only the value can be moved.
  attributes_.emplace(Utility::toupper(Utility::trimEndWhiteSpace(key)), std::move(value));

}

std::string kgl::Attributes::getHGNC() const {

  std::string hgnc_id;
  constexpr size_t HGNC_LEN = std::string_view{HGNC_}.size();

  auto const [first, last] = attributes_.equal_range(DBXREF_);
  // Keep scanning the whole DBXREF range so the LAST matching HGNC attribute wins,
  // preserving the reference semantics (no early break).
  for (auto const& [key, attrib] : std::ranges::subrange(first, last)) {

    if (attrib.starts_with(HGNC_)) {

      hgnc_id = attrib.substr(HGNC_LEN);
      hgnc_id = Utility::trimEndWhiteSpace(hgnc_id);

    }

  }

  return hgnc_id;

}

void kgl::Attributes::getSuperFeatureIds(std::vector<std::string> &value_vec) const {

  value_vec.clear();
  std::vector<std::string> super_vec;
  getAttributes(SUPER_FEATURE_KEY, super_vec);

  for (auto const& super_feature : super_vec) {

    auto parsed_features = Utility::charTokenizer(super_feature, SUPER_FEATURE_DELIMITER);
    std::ranges::copy(parsed_features, std::back_inserter(value_vec));

  }

}
