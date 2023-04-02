//
// Created by kellerberrin on 9/7/21.
//


#include "kgl_runtime_resource.h"

namespace kgl = kellerberrin::genome;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container object to hold resources as they are passed to analysis packages.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::vector<std::shared_ptr<const kgl::ResourceBase>> kgl::AnalysisResources::getResources( RuntimeResourceType resource,
                                                                                            const std::string& resource_ident) const {

  std::vector<std::shared_ptr<const ResourceBase>> resource_vector;
  for (auto const& [resource_type, resource_ptr] :  resource_map_) {

    if (not resource_ident.empty()) {

      if (resource == resource_type and resource_ptr->identifier() == resource_ident) {

        resource_vector.push_back(resource_ptr);

      }

    } else {

      if (resource == resource_type) {

        resource_vector.push_back(resource_ptr);

      }

    }

  }

  return resource_vector;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects encode XML resource definitions and pass then to the resource constructor.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


std::optional<const std::string> kgl::ResourceParameters::getParameter(const std::string& parameter_key) const {

  if (parameter_map_.contains(parameter_key)) {

    auto const& [param_key, param_value] = *parameter_map_.find(parameter_key);
    return param_value;

  }

  return std::nullopt;

}


void kgl::ResourceParameters::setParameter(const std::string& parameter_key, const std::string& parameter_value) {

  if (parameter_map_.contains(parameter_key)) {

    auto old_value = parameter_map_[parameter_key];
    ExecEnv::log().warn("ResourceParameters::setParameter; parameter key: '{}' already exists; current value: '{}', new value: '{}'",
                        parameter_key, old_value, parameter_value);

  }

  parameter_map_[parameter_key] = parameter_value;

}


std::optional<kgl::ResourceParameters> kgl::ResourceDefinitions::retrieve(const std::string& resource_type,
                                                                          const std::string& resource_ident) const {

  auto const [lower_iter, upper_iter] = resource_definition_map_.equal_range(resource_type);
  for (auto iter = lower_iter; iter != upper_iter; ++iter) {

    auto const& [type, resource_parameters] = *iter;
    if (resource_parameters.resourceIdent() == resource_ident) {

      return resource_parameters;

    }

  }

  return std::nullopt;

}
