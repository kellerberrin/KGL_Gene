//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_pf3kcoi.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_pf3k_coi.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view PF3K_COI_IDENT{"Pf3KCOIIdent"};
  constexpr std::string_view PF3K_COI_FILE{"Pf3KCOIFile"};
}


/// Returns the Pf3K COI resource type identifier.
std::string_view kgl::Pf3KCOIResourceFactory::resourceType() const {

  return ResourceType::PF3K_COI;

}


/// Parse the Pf3K COI resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::Pf3KCOIResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "Pf3KCOIResourceFactory::parseXML");
  reader.requiredString(PF3K_COI_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::PF3K_COI), ident);
  std::string file_name;
  reader.requiredFile(PF3K_COI_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(PF3K_COI_FILE), std::move(file_name));
  return params;

}


/// Construct the Pf3K COI resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::Pf3KCOIResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(PF3K_COI_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("Pf3KCOIResourceFactory::load; COI file not defined for: {}",
                            parameters.resourceIdent());

  }

  Pf3kCOIParser parser;
  if (not parser.parseCOIPf3k(file_opt.value())) {

    ExecEnv::log().critical("Pf3KCOIResourceFactory::load; failed to parse file: {}", file_opt.value());

  }

  return std::make_shared<Pf3kCOIResource>(parameters.resourceIdent(), parser.parsedText());

}
