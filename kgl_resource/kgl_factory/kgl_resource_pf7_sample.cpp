//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_pf7_sample.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_pf7_sample_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view PF7SAMPLE_IDENT{"Pf7SampleIdent"};
  constexpr std::string_view PF7SAMPLE_FILE{"Pf7SampleFile"};
}


/// Returns the Pf7 sample data resource type identifier.
std::string_view kgl::Pf7SampleResourceFactory::resourceType() const {

  return ResourceType::PF7_SAMPLE;

}


/// Parse the Pf7 sample resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::Pf7SampleResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "Pf7SampleResourceFactory::parseXML");
  reader.requiredString(PF7SAMPLE_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::PF7_SAMPLE), ident);
  std::string file_name;
  reader.requiredFile(PF7SAMPLE_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(PF7SAMPLE_FILE), std::move(file_name));
  return params;

}


/// Construct the Pf7 sample resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::Pf7SampleResourceFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(PF7SAMPLE_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("Pf7SampleResourceFactory::load; sample file not defined for: {}",
                            parameters.resourceIdent());

  }

  ParsePf7Sample parser;
  if (not parser.parsePf7SampleFile(file_opt.value())) {

    ExecEnv::log().critical("Pf7SampleResourceFactory::load; failed to parse file: {}", file_opt.value());

  }

  return std::make_shared<Pf7SampleResource>(parameters.resourceIdent(), parser.getPf7SampleVector());

}
