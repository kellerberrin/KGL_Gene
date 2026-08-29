//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_xml_parsers.h"
#include "kgl_resource_factory.h"
#include "kgl_xml_detail.h"


namespace kgl = kellerberrin::genome;
namespace kglx = kellerberrin::genome::xml;

namespace {

  constexpr std::string_view RESOURCE_LIST_{"resourceList"};

}


/// Parse resource definitions from the runtime XML configuration.
kgl::ResourceDefinitions kglx::parseResources(const PropertyTree& tree, std::string_view work_directory) {

  ResourceDefinitions resource_definitions;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not tree.getPropertyTreeVector(runtimeKey({RESOURCE_LIST_}), property_tree_vector)) {

    ExecEnv::log().warn("parseResources; No runtime resources specified.");
    return resource_definitions;

  }

  for (auto const& [tree_type, sub_tree] : property_tree_vector) {

    const ResourceFactory* factory = findResourceFactory(tree_type);
    if (not factory) {

      ExecEnv::log().warn("parseResources; Unknown resource type '{}'; ignored", tree_type);
      continue;

    }

    auto resource_opt = factory->parseXML(sub_tree, work_directory);
    if (resource_opt) {

      resource_definitions.insert({std::string(factory->resourceType()), resource_opt.value()});

    }

  }

  return resource_definitions;

}