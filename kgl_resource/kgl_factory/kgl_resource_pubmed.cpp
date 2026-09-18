//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_pubmed.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_pubmed_resource.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view PUBMED_IDENT{"pubmedApiIdent"};
  constexpr std::string_view PUBMED_PUBLICATION_CACHE{"pubmedPublicationCache"};
  constexpr std::string_view PUBMED_CITATION_CACHE{"pubmedCitationCache"};
}


/// Returns the PubMed API resource type identifier.
std::string_view kgl::PubmedResourceFactory::resourceType() const {

  return ResourceType::PUBMED_API;

}


/// Parse the PubMed API resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::PubmedResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "PubmedResourceFactory::parseXML");
  reader.requiredString(PUBMED_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::PUBMED_API), ident);
  std::string publication_file, citation_file;
  reader.createFile(PUBMED_PUBLICATION_CACHE, work_directory, publication_file)
        .createFile(PUBMED_CITATION_CACHE, work_directory, citation_file);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(PUBMED_PUBLICATION_CACHE), std::move(publication_file));
  params.setParameter(std::string(PUBMED_CITATION_CACHE), std::move(citation_file));
  return params;

}


/// Construct the PubMed API resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::PubmedResourceFactory::load(const ResourceParameters& parameters) const {

  auto publication_opt = parameters.getParameter(std::string(PUBMED_PUBLICATION_CACHE));
  auto citation_opt = parameters.getParameter(std::string(PUBMED_CITATION_CACHE));
  if (not publication_opt or not citation_opt) {

    ExecEnv::log().critical("PubmedResourceFactory::load; missing cache file parameters for: {}",
                            parameters.resourceIdent());

  }

  return std::make_shared<PubmedRequester>(parameters.resourceIdent(),
                                           publication_opt.value(),
                                           citation_opt.value());

}
