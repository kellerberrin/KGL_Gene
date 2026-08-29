//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_aux.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_hsgenome_aux.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view AUX_IDENT{"auxIdent"};
  constexpr std::string_view AUX_FILE{"auxFile"};
}


/// Returns the genome auxiliary information resource type identifier.
std::string_view kgl::AuxResourceFactory::resourceType() const {

  return ResourceType::GENOME_AUX;

}


/// Parse the genome auxiliary resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::AuxResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "AuxResourceFactory::parseXML");
  reader.requiredString(AUX_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::GENOME_AUX), ident);
  std::string file_name;
  reader.requiredFile(AUX_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(AUX_FILE), std::move(file_name));
  return params;

}


/// Construct the genome auxiliary resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::AuxResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(AUX_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("AuxResourceFactory::load; aux file not defined for: {}",
                            parameters.resourceIdent());

  }

  auto genome_aux_data = std::make_shared<HsGenomeAuxData>(parameters.resourceIdent());
  ParseHsGenomeAuxFile genome_aux_parser(genome_aux_data);
  genome_aux_parser.readParseImpl(file_opt.value());

  return genome_aux_data;

}
