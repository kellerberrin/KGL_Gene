//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_ontology.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_ontology_database.h"

namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

namespace {
  constexpr std::string_view ONTOLOGY_IDENT{"ontologyIdent"};
  constexpr std::string_view GAF_ANNOTATION_FILE{"gafFile"};
  constexpr std::string_view GO_ONTOLOGY_FILE{"goFile"};
}


/// Returns the ontology database resource type identifier.
std::string_view kgl::OntologyResourceFactory::resourceType() const {

  return ResourceType::ONTOLOGY;

}


/// Parse the ontology resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::OntologyResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "OntologyResourceFactory::parseXML");
  reader.requiredString(ONTOLOGY_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::ONTOLOGY), ident);
  std::string gaf, go;
  reader.requiredFile(GAF_ANNOTATION_FILE, work_directory, gaf)
        .requiredFile(GO_ONTOLOGY_FILE, work_directory, go);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(GAF_ANNOTATION_FILE), std::move(gaf));
  params.setParameter(std::string(GO_ONTOLOGY_FILE), std::move(go));
  return params;

}


/// Construct the ontology database resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::OntologyResourceFactory::load(const ResourceParameters& parameters) const {

  auto gaf_opt = parameters.getParameter(std::string(GAF_ANNOTATION_FILE));
  auto go_opt = parameters.getParameter(std::string(GO_ONTOLOGY_FILE));
  if (not gaf_opt or not go_opt) {

    ExecEnv::log().critical("OntologyResourceFactory::load; missing required parameters for: {}",
                            parameters.resourceIdent());

  }

  return std::make_shared<const kol::OntologyDatabase>(parameters.resourceIdent(),
                                                        go_opt.value(),
                                                        gaf_opt.value());

}