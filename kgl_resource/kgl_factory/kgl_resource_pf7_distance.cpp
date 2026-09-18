//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_pf7_distance.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_pf7_genetic_distance_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view PF7DISTANCE_IDENT{"Pf7DistanceIdent"};
  constexpr std::string_view PF7DISTANCE_MATRIXFILE{"Pf7DistanceMatrixFile"};
  constexpr std::string_view PF7DISTANCE_IDFILE{"Pf7DistanceIDFile"};
}


/// Returns the Pf7 genetic distance resource type identifier.
std::string_view kgl::Pf7DistanceResourceFactory::resourceType() const {

  return ResourceType::PF7_DISTANCE;

}


/// Parse the Pf7 genetic distance resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::Pf7DistanceResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "Pf7DistanceResourceFactory::parseXML");
  reader.requiredString(PF7DISTANCE_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::PF7_DISTANCE), ident);
  std::string matrix_file, id_file;
  reader.requiredFile(PF7DISTANCE_MATRIXFILE, work_directory, matrix_file)
        .requiredFile(PF7DISTANCE_IDFILE, work_directory, id_file);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(PF7DISTANCE_MATRIXFILE), std::move(matrix_file));
  params.setParameter(std::string(PF7DISTANCE_IDFILE), std::move(id_file));
  return params;

}


/// Construct the Pf7 genetic distance resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::Pf7DistanceResourceFactory::load(const ResourceParameters& parameters) const {

  auto matrix_opt = parameters.getParameter(std::string(PF7DISTANCE_MATRIXFILE));
  auto id_opt = parameters.getParameter(std::string(PF7DISTANCE_IDFILE));
  if (not matrix_opt or not id_opt) {

    ExecEnv::log().critical("Pf7DistanceResourceFactory::load; missing required parameters for: {}",
                            parameters.resourceIdent());

  }

  ParsePf7Distance parser;
  if (not parser.parsePf7Distance(matrix_opt.value(), id_opt.value())) {

    ExecEnv::log().critical("Pf7DistanceResourceFactory::load; failed to parse matrix: {}, id file: {}",
                            matrix_opt.value(), id_opt.value());

  }

  return std::make_shared<Pf7GeneticDistanceResource>(parameters.resourceIdent(),
                                                      parser.getSampleMap(),
                                                      parser.getDistanceMatrix());

}
