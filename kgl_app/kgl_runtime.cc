//
// Created by kellerberrin on 29/4/20.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_runtime.h"


namespace kgl = kellerberrin::genome;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

kgl::DataFileParserEnum kgl::BaseFileInfo::getParserType(const std::string& parser_type) const {

  std::string parser_upper = Utility::toupper(parser_type);

  for (auto const& [parser_type, parser_string] : implementated_parsers_) {

    if (parser_upper == Utility::toupper(parser_string)) {

      return parser_type;

    }

  }

  return DataFileParserEnum::NotImplemented;

}


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
