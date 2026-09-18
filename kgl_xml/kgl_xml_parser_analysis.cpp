//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_xml_reader.h"
#include "kgl_xml_detail.h"

#include <string>
#include <vector>


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  constexpr std::string_view ANALYSIS_LIST_{"analysisList"};
  constexpr std::string_view ANALYSIS_{"analysis"};
  constexpr std::string_view ANALYSIS_IDENT_{"analysisIdent"};
  constexpr std::string_view PARAMETER_RUNTIME_{"parameterRuntime"};
  constexpr std::string_view PARAMETER_BLOCK_{"parameterBlock"};

}


/// Parse the analysis map from the runtime XML configuration.
kgl::RuntimeAnalysisMap kglx::parseAnalysisMap(const PropertyTree& tree) {

  RuntimeAnalysisMap analysis_map;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({ANALYSIS_LIST_}), property_tree_vector)) {

    ExecEnv::log().info("parseAnalysisMap; no Analysis specified (NULL analysis available)");
    return analysis_map;

  }

  for (auto const& [tag, sub_tree] : property_tree_vector) {

    if (tag != ANALYSIS_) continue;

    XmlReader reader(sub_tree, "parseAnalysisMap");

    std::string analysis_ident;
    reader.requiredString(ANALYSIS_IDENT_, analysis_ident);
    if (not reader.ok()) {

      continue;

    }

    std::string param_key = std::string(PARAMETER_RUNTIME_) + "." + std::string(ACTIVE_);
    std::vector<SubPropertyTree> parameter_tree_vector;
    if (not sub_tree.getPropertyTreeVector(param_key, parameter_tree_vector)) {

      ExecEnv::log().info("parseAnalysisMap; no Parameters Specified for Analysis: {}", analysis_ident);

    }

    RuntimeParameterMap parameter_map;
    for (auto const& [param_tag, param_sub_tree] : parameter_tree_vector) {

      if (param_tag == PARAMETER_BLOCK_) {

        std::string parameter_ident = param_sub_tree.getValue();
        if (parameter_ident.empty()) {

          ExecEnv::log().error("parseAnalysisMap; No Parameter Identifier Found for Analysis: {}", analysis_ident);
          continue;

        }

        parameter_map.push_back(parameter_ident);

      }

    }

    auto [iter, result] = analysis_map.emplace(analysis_ident, RuntimeAnalysis(analysis_ident, parameter_map));
    if (not result) {

      ExecEnv::log().error("parseAnalysisMap; Could not add Analysis Ident: {} to map (duplicate)", analysis_ident);

    }

  }

  return analysis_map;

}