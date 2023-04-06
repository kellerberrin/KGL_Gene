//
// Created by kellerberrin on 5/04/23.
//


#include "kgl_properties_resource.h"
#include "kel_utility.h"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// P. Falciparum resources
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kgl = kellerberrin::genome;


std::optional<kgl::ResourceParameters> kgl::ResourceProperties::Pf7SampleResourceXML(const PropertyTree& sub_tree) const {

  std::string Pf7Sample_ident;
  if (not sub_tree.getProperty(PF7SAMPLE_IDENT_, Pf7Sample_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(PF7SAMPLE_RESOURCE_ID_, Pf7Sample_ident);

  std::string Pf7Sample_file_name;
  if (not sub_tree.getFileProperty(PF7SAMPLE_FILE_, workDirectory(), Pf7Sample_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample data file name information, ident: {}", Pf7Sample_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PF7SAMPLE_FILE_, Pf7Sample_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::ResourceProperties::Pf7FwsResourceXML(const PropertyTree& sub_tree) const {

  std::string Pf7Fws_ident;
  if (not sub_tree.getProperty(PF7FWS_IDENT_, Pf7Fws_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 FWS Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(PF7FWS_RESOURCE_ID_, Pf7Fws_ident);

  std::string Pf7Fws_file_name;
  if (not sub_tree.getFileProperty(PF7FWS_FILE_, workDirectory(), Pf7Fws_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 FWS file name information, ident: {}", Pf7Fws_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PF7FWS_FILE_, Pf7Fws_file_name);

  return resource_parameters;

}


std::optional<kgl::ResourceParameters> kgl::ResourceProperties::Pf7DistanceResourceXML(const PropertyTree& sub_tree) const {

  std::string Pf7Distance_ident;
  if (not sub_tree.getProperty(PF7DISTANCE_IDENT_, Pf7Distance_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 FWS Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(PF7DISTANCE_RESOURCE_ID_, Pf7Distance_ident);

  std::string Pf7Matrix_file_name;
  if (not sub_tree.getFileProperty(PF7DISTANCE_MATRIXFILE_, workDirectory(), Pf7Matrix_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 Distance Matrix file name information, ident: {}", Pf7Distance_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PF7DISTANCE_MATRIXFILE_, Pf7Matrix_file_name);

  std::string Pf7ID_file_name;
  if (not sub_tree.getFileProperty(PF7DISTANCE_IDFILE_, workDirectory(), Pf7ID_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 Distance sample ID file name information, ident: {}", Pf7Distance_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PF7DISTANCE_IDFILE_, Pf7ID_file_name);

  return resource_parameters;

}

