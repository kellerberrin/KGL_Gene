//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_entrez.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_entrez_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view ENTREZ_IDENT{"entrezIdent"};
  constexpr std::string_view ENTREZ_FILE{"entrezFile"};
}


/// Returns the Entrez gene resource type identifier.
std::string_view kgl::EntrezResourceFactory::resourceType() const {

  return ResourceType::ENTREZ;

}


/// Parse the Entrez gene resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::EntrezResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "EntrezResourceFactory::parseXML");
  reader.requiredString(ENTREZ_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::ENTREZ), ident);
  std::string file_name;
  reader.requiredFile(ENTREZ_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(ENTREZ_FILE), std::move(file_name));
  return params;

}


/// Construct the Entrez gene resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::EntrezResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(ENTREZ_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("EntrezResourceFactory::load; entrez file not defined for: {}",
                            parameters.resourceIdent());

  }

  ParseEntrez entrez_parser;
  if (not entrez_parser.parseEntrezFile(file_opt.value())) {

    ExecEnv::log().critical("EntrezResourceFactory::load; failed to parse file: {}", file_opt.value());

  }

  return std::make_shared<EntrezResource>(parameters.resourceIdent(), entrez_parser.getEntrezVector());

}
