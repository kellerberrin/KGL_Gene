//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_citation.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_citation_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view CITATION_IDENT{"citationIdent"};
  constexpr std::string_view CITATION_FILE{"citationFile"};
}


/// Returns the allele citation resource type identifier.
std::string_view kgl::CitationResourceFactory::resourceType() const {

  return ResourceType::CITATION;

}


/// Parse the citation resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::CitationResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "CitationResourceFactory::parseXML");
  reader.requiredString(CITATION_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::CITATION), ident);
  std::string file_name;
  reader.requiredFile(CITATION_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(CITATION_FILE), std::move(file_name));
  return params;

}


/// Construct the allele citation resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::CitationResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(CITATION_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("CitationResourceFactory::load; citation file not defined for: {}",
                            parameters.resourceIdent());

  }

  ParseCitations citation_parser;
  if (not citation_parser.parseCitationFile(file_opt.value())) {

    ExecEnv::log().critical("CitationResourceFactory::load; failed to parse file: {}", file_opt.value());

  }

  return std::make_shared<CitationResource>(parameters.resourceIdent(), citation_parser.getCitationMap());

}
