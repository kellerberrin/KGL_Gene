//
// Created by kellerberrin on 29/4/20.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_runtime.h"


namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

const kgl::ContigId_t& kgl::ContigAliasMap::lookupAlias(const ContigId_t& alias) const {

  auto result = alias_map_.find(alias);

  if (result == alias_map_.end()) {

    ExecEnv::log().error("ContigAliasMap::lookupAlias(); Alias: {} Not Found", alias);
    return alias;

  }

  return result->second;

}

void kgl::ContigAliasMap::setAlias(const ContigId_t& alias, const ContigId_t& contig_id) {

  auto result = alias_map_.insert(std::pair<std::string, std::string>(alias, contig_id));

  if (not result.second) {

    ExecEnv::log().error("ContigAliasMap::setAlias(); Cannot register Alias: {} for Contig: {} (duplicate)", alias, contig_id);

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

std::optional<const kgl::EvidenceInfoSet> kgl::VariantEvidenceMap::lookupEvidence(const std::string& evidence_ident) const {

  auto result = evidence_map_.find(evidence_ident);

  if (result == evidence_map_.end()) {

    ExecEnv::log().error("VariantEvidenceMap::lookupEvidence(); Cannot find evidence ident: {}", evidence_ident);
    return std::nullopt;

  }

  return result->second;

}


void kgl::VariantEvidenceMap::setEvidence(const std::string& evidence_ident, const std::set<std::string>& info_list) {

  auto result = evidence_map_.emplace(evidence_ident, info_list);

  if (not result.second) {

    ExecEnv::log().error("VariantEvidenceMap::setEvidence(); Cannot register evidence ident: {} (duplicate)", evidence_ident);

  }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

bool kgl::ActiveParameterList::addNamedParameterVector(const NamedParameterVector& named_vector) {

  auto insert_result = active_parameter_vectors_.try_emplace(named_vector.first, named_vector);
  auto [iterator, result] = insert_result;

  if (not result) {

    ExecEnv::log().error("ActiveParameterList::addNamedParameterVector; could not add named parameter vector: {} (duplicate)", named_vector.first);

  }

  return result;

}


kgl::ActiveParameterList kgl::ActiveParameterList::createParameterList(const std::vector<std::string>& active_idents) const {

  ActiveParameterList active_parameters;

  for (auto const& ident : active_idents) {

    auto result = active_parameter_vectors_.find(ident);
    if (result == active_parameter_vectors_.end()) {

      ExecEnv::log().error("ActiveParameterList::createParameterList; specified active parameter: {} not found in master list", ident);

    } else {

      auto const& [ident, named_param_vector] = *result;
      active_parameters.addNamedParameterVector(named_param_vector);

    }

  }

  return active_parameters;

}


std::vector<std::string> kgl::ParameterMap::retrieve(const std::string& ident) const {

  std::vector<std::string> values;

  auto range_iterator = parameter_map_.equal_range(ident);

  for (auto it = range_iterator.first; it != range_iterator.second; ++it)
  {

    auto [key, value] = *it;
    values.push_back(value);

  }

  return values;

}


std::optional<std::vector<double>> kgl::ParameterMap::getFloat(const std::string& ident, size_t vec_size) const {

  std::vector<double> value_vector;

  std::vector<std::string> string_vector = retrieve(ident);

  if (string_vector.size() != vec_size and vec_size != ANY_SIZE) {

    return std::nullopt;

  }

  for (auto const& str_value : string_vector) {

    try {

      value_vector.push_back(std::stod(str_value));

    } catch(...) {

      ExecEnv::log().error("ParameterMap::getFloat; parameter ident: {}, has invalid double value: {}", ident, str_value);
      return std::nullopt;

    }

  }


  return value_vector;

}


std::optional<std::vector<std::string>> kgl::ParameterMap::getString(const std::string& ident, size_t vec_size) const {

  std::vector<std::string> value_vector = retrieve(ident);

  if (value_vector.size() != vec_size and vec_size != ANY_SIZE) {

    return std::nullopt;

  }

  return value_vector;

}


std::optional<std::vector<int64_t>> kgl::ParameterMap::getInteger(const std::string& ident, size_t vec_size) const {

  std::vector<int64_t> value_vector;

  std::vector<std::string> string_vector = retrieve(ident);

  if (string_vector.size() != vec_size and vec_size != ANY_SIZE) {

    return std::nullopt;

  }

  for (auto const& str_value : string_vector) {

    try {

      value_vector.push_back(std::stoll(str_value));

    } catch(...) {

      ExecEnv::log().error("ParameterMap::getInteger; parameter ident: {}, has invalid signed integer value: {}", ident, str_value);
      return std::nullopt;

    }

  }

  return value_vector;

}



std::optional<std::vector<size_t>> kgl::ParameterMap::getSize(const std::string& ident, size_t vec_size) const {

  std::vector<size_t> value_vector;

  std::vector<std::string> string_vector = retrieve(ident);

  if (string_vector.size() != vec_size and vec_size != ANY_SIZE) {

    return std::nullopt;

  }

  for (auto const& str_value : string_vector) {

    try {

      value_vector.push_back(std::stoull(str_value));

    } catch(...) {

      ExecEnv::log().error("ParameterMap::getSize; parameter ident: {}, has invalid unsigned integer value: {}", ident, str_value);
      return std::nullopt;

    }

  }

  return value_vector;

}


std::optional<bool> kgl::ParameterMap::getBool(const std::string& ident) const {

  static const std::string true_str{"TRUE"};
  static const std::string false_str{"FALSE"};

  std::vector<std::string> string_vector = retrieve(ident);

  if (string_vector.size() != 1) {

    return std::nullopt;

  }

  auto bool_str = Utility::toupper(string_vector.front());

  if (bool_str == true_str) {

    return true;

  } else if (bool_str == false_str) {

    return false;

  }

  ExecEnv::log().error("ParameterMap::getBool; parameter ident: {}, has invalid boolean value: {}", ident, string_vector.front());

  return std::nullopt;

}
