//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_pf7_fws.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_pf7_fws_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view PF7FWS_IDENT{"Pf7FwsIdent"};
  constexpr std::string_view PF7FWS_FILE{"Pf7FwsFile"};
}


/// Returns the Pf7 FWS resource type identifier.
std::string_view kgl::Pf7FwsResourceFactory::resourceType() const {

  return ResourceType::PF7_FWS;

}


/// Parse the Pf7 FWS resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::Pf7FwsResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "Pf7FwsResourceFactory::parseXML");
  reader.requiredString(PF7FWS_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::PF7_FWS), ident);
  std::string file_name;
  reader.requiredFile(PF7FWS_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(PF7FWS_FILE), std::move(file_name));
  return params;

}


/// Construct the Pf7 FWS resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::Pf7FwsResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(PF7FWS_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("Pf7FwsResourceFactory::load; FWS file not defined for: {}",
                            parameters.resourceIdent());

  }

  ParsePf7Fws parser;
  if (not parser.parsePf7FwsFile(file_opt.value())) {

    ExecEnv::log().critical("Pf7FwsResourceFactory::load; failed to parse file: {}", file_opt.value());

  }

  return std::make_shared<Pf7FwsResource>(parameters.resourceIdent(), parser.getPf7FwsVector());

}
