//
// Created by kellerberrin on 9/7/21.
//


#include "kgl_resource_db.h"

namespace kgl = kellerberrin::genome;





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


std::string  kgl::AnalysisResources::resourceDescription(RuntimeResourceType resource) {

  for (auto const& [resource_key, description] : resource_description) {

    if (resource_key == resource) {

      return description;

    }

  }

  return "Resource Description Not Found";

}


