//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_xml_reader.h"
#include "kgl_xml_detail.h"

#include <string>
#include <utility>
#include <vector>


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  constexpr std::string_view PACKAGE_LIST_{"packageList"};
  constexpr std::string_view PACKAGE_{"package"};
  constexpr std::string_view PACKAGE_IDENT_{"packageIdent"};
  constexpr std::string_view PACKAGE_ANALYSIS_LIST_{"analysisList"};
  constexpr std::string_view ANALYSIS_{"analysis"};
  constexpr std::string_view PACKAGE_RESOURCE_LIST_{"resourceList"};
  constexpr std::string_view PACKAGE_ITERATION_{"iteration"};
  constexpr std::string_view PACKAGE_ITERATION_LIST_{"iterationList"};
  constexpr std::string_view DATA_FILE_IDENT_{"dataFileIdent"};

}


/// Parse the package map from the runtime XML configuration.
kgl::RuntimePackageMap kglx::parsePackageMap(const PropertyTree& tree) {

  RuntimePackageMap package_map;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({PACKAGE_LIST_}), property_tree_vector)) {

    ExecEnv::log().error("parsePackageMap; no Analysis Packages Specified (must be at least one)");
    return package_map;

  }

  for (auto const& [tag, sub_tree] : property_tree_vector) {

    if (tag != PACKAGE_) continue;

    XmlReader reader(sub_tree, "parsePackageMap");

    std::string package_ident;
    reader.requiredString(PACKAGE_IDENT_, package_ident);
    if (not reader.ok()) {

      continue;

    }

    // Analysis list.
    std::vector<SubPropertyTree> analysis_sub_trees;
    if (not sub_tree.getPropertyTreeVector(std::string(PACKAGE_ANALYSIS_LIST_), analysis_sub_trees)) {

      ExecEnv::log().warn("parsePackageMap; No Analysis specified for Package: {}", package_ident);

    }
    std::vector<std::string> analysis_vector;
    for (auto const& [analysis_tag, analysis_sub_tree] : analysis_sub_trees) {

      if (analysis_tag == ANALYSIS_) {

        std::string analysis_ident = analysis_sub_tree.getValue();
        if (analysis_ident.empty()) {

          ExecEnv::log().critical("parsePackageMap; Package: {}, No Analysis identifier specified", package_ident);

        } else {

          analysis_vector.push_back(analysis_ident);

        }

      }

    }

    // Resource list.
    std::vector<SubPropertyTree> resource_vector;
    if (not sub_tree.getPropertyTreeVector(std::string(PACKAGE_RESOURCE_LIST_), resource_vector)) {

      ExecEnv::log().warn("parsePackageMap; No resources specified for Package: {}", package_ident);

    }
    std::vector<std::pair<std::string, std::string>> resources_def;
    for (auto const& [resource_type, resource_sub_tree] : resource_vector) {

      std::string resource_identifier = resource_sub_tree.getValue();
      if (resource_identifier.empty()) {

        ExecEnv::log().critical("parsePackageMap; Package: {}, Resource: {} No resource identifier specified",
                                package_ident, resource_type);

      } else {

        resources_def.emplace_back(resource_type, resource_identifier);

      }

    }

    // Iteration list (VCF files).
    std::vector<SubPropertyTree> iteration_vector;
    if (not sub_tree.getPropertyTreeVector(std::string(PACKAGE_ITERATION_LIST_), iteration_vector)) {

      ExecEnv::log().warn("parsePackageMap; No Iteration items (VCF Files) specified for Package: {}", package_ident);

    }
    std::vector<std::vector<std::string>> vector_iteration_files;
    for (auto const& [iteration_tag, iteration_sub_tree] : iteration_vector) {

      if (iteration_tag != PACKAGE_ITERATION_) continue;

      std::vector<std::string> iteration_files;
      if (not iteration_sub_tree.getNodeVector(std::string(DATA_FILE_IDENT_), iteration_files)) {

        ExecEnv::log().critical("parsePackageMap; No Iteration items (VCF Files) specified for Package: {}", package_ident);

      }

      if (not iteration_files.empty()) {

        vector_iteration_files.push_back(iteration_files);

      }

    }

    auto [iter, result] = package_map.emplace(package_ident,
                                               RuntimePackage(package_ident, analysis_vector,
                                                              resources_def, vector_iteration_files));
    if (not result) {

      ExecEnv::log().error("parsePackageMap; Could not add Package Ident: {} to map (duplicate)", package_ident);

    }

  }

  return package_map;

}