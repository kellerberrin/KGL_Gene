//
// Created by kellerberrin on 29/4/20.
//

#include "kel_exec_env.h"
#include "kel_utility.h"
#include "kgl_runtime.h"


namespace kgl = kellerberrin::genome;


kgl::VCFParserEnum kgl::VCFFileInfo::getParserType(const std::string& parser_type) const {

  std::string parser_upper = Utility::toupper(parser_type);

  for (auto const& [parser_type, parser_string] : implementated_parsers_) {

    if (parser_upper == Utility::toupper(parser_string)) {

      return parser_type;

    }

  }

  return VCFParserEnum::NotImplemented;

}



const kgl::ContigId_t& kgl::ContigAliasMap::lookupAlias(const ContigId_t& alias) {

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

