//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_RESOURCE_DB_H
#define KGL_RESOURCE_DB_H


#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ResourceBase inherited by all resources.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class RuntimeResourceType { GENOME_DATABASE, ONTOLOGY_DATABASE, GENE_NOMENCLATURE };

class ResourceBase {

public:

  ResourceBase() = default;
  virtual ~ResourceBase() = default;

  [[nodiscard]] virtual RuntimeResourceType getResourceType() const = 0;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container object to hold resources as they are passed to analysis packages.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ResourceMap = std::multimap<RuntimeResourceType, std::shared_ptr<const ResourceBase>>;


class AnalysisResources {

public:

  AnalysisResources() = default;
  ~AnalysisResources() = default;

  void addResource(const std::shared_ptr<const ResourceBase>& resource_ptr) {

    resource_map_.emplace(resource_ptr->getResourceType(), resource_ptr);

  }

  [[nodiscard]] std::vector<std::shared_ptr<const ResourceBase>> getResources(RuntimeResourceType resource) const {

    std::vector<std::shared_ptr<const ResourceBase>> resource_vector;
    for (auto const& [resource_type, resource_ptr] :  resource_map_) {

      if (resource == resource_type) {

        resource_vector.push_back(resource_ptr);

      }

    }

    return resource_vector;

  }

  [[nodiscard]] const ResourceMap& getMap() const { return resource_map_; }

private:

  ResourceMap resource_map_;

};



}   // end namespace

#endif //KGL_RESOURCE_DB_H
