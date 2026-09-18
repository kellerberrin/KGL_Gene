//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_genealogy.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_hsgenealogy_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view GENEALOGY_IDENT{"genealogyIdent"};
  constexpr std::string_view GENEALOGY_FILE{"genealogyFile"};
}


/// Returns the genome genealogy resource type identifier.
std::string_view kgl::GenealogyResourceFactory::resourceType() const {

  return ResourceType::GENEALOGY;

}


/// Parse the genealogy resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::GenealogyResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "GenealogyResourceFactory::parseXML");
  reader.requiredString(GENEALOGY_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::GENEALOGY), ident);
  std::string file_name;
  reader.requiredFile(GENEALOGY_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(GENEALOGY_FILE), std::move(file_name));
  return params;

}


/// Construct the genome genealogy resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::GenealogyResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(GENEALOGY_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("GenealogyResourceFactory::load; genealogy file not defined for: {}",
                            parameters.resourceIdent());

  }

  auto genealogy_data = std::make_shared<HsGenomeGenealogyData>(parameters.resourceIdent());
  ParseHsGenomeGenealogyFile genealogy_parser(genealogy_data);
  genealogy_parser.readParseImpl(file_opt.value());

  return genealogy_data;

}
