//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_xml_reader.h"
#include "kgl_xml_detail.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  constexpr std::string_view PARAMETER_LIST_{"parameterList"};
  constexpr std::string_view PARAMETER_BLOCK_{"parameterBlock"};
  constexpr std::string_view PARAMETER_NAME_{"parameterName"};
  constexpr std::string_view PARAMETER_VECTOR_{"parameterVector"};
  constexpr std::string_view PARAMETER_{"parameter"};
  constexpr std::string_view PARAMETER_IDENT_{"parameterIdent"};
  constexpr std::string_view PARAMETER_VALUE_{"parameterValue"};

}


/// Parse named parameter blocks from the runtime XML configuration.
kgl::ActiveParameterList kglx::parseParameters(const PropertyTree& tree) {

  ActiveParameterList defined_named_parameters;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({PARAMETER_LIST_}), property_tree_vector)) {

    ExecEnv::log().warn("parseParameters; No Parameter Blocks specified");
    return defined_named_parameters;

  }

  for (auto const& [sub_tree_tag, sub_tree] : property_tree_vector) {

    if (sub_tree_tag != PARAMETER_BLOCK_) continue;

    NamedParameterVector named_parameter_vector;

    XmlReader reader(sub_tree, "parseParameters");

    std::string block_name;
    reader.requiredString(PARAMETER_NAME_, block_name);
    if (not reader.ok()) {

      ExecEnv::log().error("parseParameters; No block <parameterName> specified, block skipped");
      continue;

    }
    named_parameter_vector.first = Utility::trimEndWhiteSpace(block_name);

    std::vector<SubPropertyTree> property_block_vector;
    if (not sub_tree.getPropertySubTreeVector(property_block_vector)) {

      ExecEnv::log().warn("parseParameters; No Parameter Vectors specified");

    }

    ParameterVector parameter_vector;
    for (auto const& [block_tree_tag, block_sub_tree] : property_block_vector) {

      if (block_tree_tag != PARAMETER_VECTOR_) continue;

      std::vector<SubPropertyTree> property_vector;
      if (not block_sub_tree.getPropertySubTreeVector(property_vector)) {

        continue;

      }

      ParameterMap parameter_map;
      for (auto const& [item_tag, item_tree] : property_vector) {

        if (item_tag != PARAMETER_) continue;

        XmlReader param_reader(item_tree, "parseParameters(item)");

        std::string parameter_ident;
        std::string parameter_value;
        param_reader.requiredString(PARAMETER_IDENT_, parameter_ident)
                    .requiredString(PARAMETER_VALUE_, parameter_value);
        if (not param_reader.ok()) {

          continue;

        }

        parameter_map.insert(Utility::trimEndWhiteSpace(parameter_ident),
                             Utility::trimEndWhiteSpace(parameter_value));

      }

      parameter_vector.push_back(parameter_map);

    }

    named_parameter_vector.second = parameter_vector;
    defined_named_parameters.addNamedParameterVector(named_parameter_vector);

  }

  return defined_named_parameters;

}