//
// Created by kellerberrin on 29/4/20.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_runtime.h"


namespace kgl = kellerberrin::genome;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

/// Looks up a contig alias and returns the canonical contig identifier.
const kgl::ContigId_t& kgl::ContigAliasMap::lookupAlias(const ContigId_t& alias) const {

  auto result = alias_map_.find(alias);
  if (result == alias_map_.end()) {

    ExecEnv::log().error("ContigAliasMap::lookupAlias(); Alias: {} Not Found", alias);
    return alias;

  }

  // Unpack the alias map.
  auto const& [contig_key, alias_pair] = *result;
  auto const& [contig_alias, contig_type] = alias_pair;

  return contig_alias;

}

/// Looks up the chromosome type associated with a contig alias.
kgl::ChromosomeType kgl::ContigAliasMap::lookupType(const ContigId_t& alias) const {

  auto result = alias_map_.find(alias);
  if (result == alias_map_.end()) {

    ExecEnv::log().error("ContigAliasMap::lookupType(); Alias: {} Not Found", alias);
    return ChromosomeType::AUTOSOMAL;

  }

  // Unpack the alias map.
  auto const& [contig_key, alias_pair] = *result;
  auto const& [contig_alias, contig_type] = alias_pair;

  return contig_type;

}

/// Sets a contig alias with the specified canonical contig identifier and chromosome type string.
void kgl::ContigAliasMap::setAlias(const ContigId_t &alias, const ContigId_t &contig_id, const std::string &chromosome_type) {

  ChromosomeType chrom_type{ChromosomeType::AUTOSOMAL};
  if (chromosome_type == AUTOSOME_) {

    chrom_type = ChromosomeType::AUTOSOMAL;

  } else if (chromosome_type == ALLOSOME_X_) {

    chrom_type = ChromosomeType::ALLOSOME_X;

  } else if (chromosome_type == ALLOSOME_Y_) {

    chrom_type = ChromosomeType::ALLOSOME_Y;

  } else if (chromosome_type == MITOCHONDRIA_) {

    chrom_type = ChromosomeType::MITOCHONDRIA;

  } else {

    ExecEnv::log().error("ContigAliasMap::setAlias; unknown chromosome type: {} must one of {}, {}, {} or {}, see alias xml file.",
                          chromosome_type, AUTOSOME_, ALLOSOME_X_, ALLOSOME_Y_, MITOCHONDRIA_);

  }

  auto result = alias_map_.try_emplace(alias, std::pair<std::string, ChromosomeType>{contig_id, chrom_type});
  if (not result.second) {

    ExecEnv::log().error("ContigAliasMap::setAlias(); Cannot register Alias: {} for Contig: {} (duplicate)", alias, contig_id);

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

/// Looks up the evidence info set associated with the given evidence identifier.
std::optional<const kgl::EvidenceInfoSet> kgl::VariantEvidenceMap::lookupEvidence(const std::string& evidence_ident) const {

  auto result = evidence_map_.find(evidence_ident);

  if (result == evidence_map_.end()) {

    ExecEnv::log().error("VariantEvidenceMap::lookupEvidence(); Cannot find evidence ident: {}", evidence_ident);
    return std::nullopt;

  }

  return result->second;

}


/// Registers an evidence identifier with its associated set of info fields.
void kgl::VariantEvidenceMap::setEvidence(const std::string& evidence_ident, const std::set<std::string>& info_list) {

  auto result = evidence_map_.emplace(evidence_ident, info_list);

  if (not result.second) {

    ExecEnv::log().error("VariantEvidenceMap::setEvidence(); Cannot register evidence ident: {} (duplicate)", evidence_ident);

  }

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

/// Adds a named parameter vector to the active parameter list.
bool kgl::ActiveParameterList::addNamedParameterVector(const NamedParameterVector& named_vector) {

  auto insert_result = active_parameter_vectors_.try_emplace(named_vector.first, named_vector);
  auto [iterator, result] = insert_result;

  if (not result) {

    ExecEnv::log().error("ActiveParameterList::addNamedParameterVector; could not add named parameter vector: {} (duplicate)", named_vector.first);

  }

  return result;

}


/// Creates a filtered parameter list containing only the named parameter vectors matching the specified identifiers.
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


/// Retrieves all values associated with the specified parameter identifier from the parameter map.
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


/// Retrieves a vector of floating-point values for the specified parameter identifier.
std::optional<std::vector<double>> kgl::ParameterMap::getFloat(const std::string& ident, size_t vec_size) const {
  return convertParameter<double>(ident, vec_size, [](const auto& s) { return std::stod(s); });
}


/// Retrieves a vector of string values for the specified parameter identifier.
std::optional<std::vector<std::string>> kgl::ParameterMap::getString(const std::string& ident, size_t vec_size) const {

  std::vector<std::string> value_vector = retrieve(ident);

  if (value_vector.size() != vec_size and vec_size != ANY_SIZE) {

    return std::nullopt;

  }

  return value_vector;

}


/// Retrieves a vector of 64-bit integer values for the specified parameter identifier.
std::optional<std::vector<int64_t>> kgl::ParameterMap::getInteger(const std::string& ident, size_t vec_size) const {
  return convertParameter<int64_t>(ident, vec_size, [](const auto& s) { return std::stoll(s); });
}



/// Retrieves a vector of size_t values for the specified parameter identifier.
std::optional<std::vector<size_t>> kgl::ParameterMap::getSize(const std::string& ident, size_t vec_size) const {
  return convertParameter<size_t>(ident, vec_size, [](const auto& s) { return std::stoull(s); });
}


/// Retrieves a boolean value for the specified parameter identifier.
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