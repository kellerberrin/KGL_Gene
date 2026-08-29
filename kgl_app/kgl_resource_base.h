//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_RESOURCE_BASE_H
#define KGL_RESOURCE_BASE_H

#include "kel_exec_env.h"

#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Base class inherited by all resources.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// Abstract base class for all loaded resources. Provides type and identifier metadata.
class ResourceBase {

public:

  ResourceBase(std::string resource_type, std::string resource_ident)
      : resource_type_(std::move(resource_type)), resource_ident_(std::move(resource_ident)) {}
  virtual ~ResourceBase() = default;

  /// Returns the resource type string.
  [[nodiscard]] const std::string& resourceType() const { return resource_type_; }
  /// Returns the resource identifier string.
  [[nodiscard]] const std::string& resourceIdent() const { return resource_ident_; }

private:

  const std::string resource_type_;
  const std::string resource_ident_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Simple container object to hold resources as they are passed to analysis packages.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ResourceMap = std::multimap<std::string, std::shared_ptr<const ResourceBase>>;

/// Container object that collects loaded resources and provides typed lookup for analysis objects.
class AnalysisResources {

public:

  AnalysisResources() = default;
  ~AnalysisResources() = default;

  /// Adds a resource to the container, keyed by resource type.
  void addResource(const std::shared_ptr<const ResourceBase>& resource_ptr) {

    resource_map_.emplace(resource_ptr->resourceType(), resource_ptr);

  }



  /// Returns all resources matching the given type (and optionally identifier).
  [[nodiscard]] std::vector<std::shared_ptr<const ResourceBase>> getResources( const std::string& resource_type,
                                                                               const std::string& resource_ident = "") const;

  /// Returns the underlying multimap of all resources.
  [[nodiscard]] const ResourceMap& getMap() const { return resource_map_; }


  /// Returns exactly one resource, cast to the specified type. Calls critical() on failure.
  template <class ResourceClass>
  [[nodiscard]] std::shared_ptr<const ResourceClass> getSingleResource(std::string resource_type, std::string resource_ident = "") const {

    auto resource_vector = getResources(resource_type, resource_ident);
    if (resource_vector.size() != 1) {

      ExecEnv::log().critical( "Request Resource Type: {} Ident: {} expected 1 resource, found: {} resources - unrecoverable error",
                               resource_type, resource_ident, resource_vector.size());

    }

    auto resource_ptr = std::dynamic_pointer_cast<const ResourceClass>(resource_vector.front());

    if (not resource_ptr) {

      ExecEnv::log().critical( "Request Resource Ident: {} invalid resource type found - unrecoverable error", resource_ident);

    }

    return resource_ptr;

  }


private:

  ResourceMap resource_map_;

};


}   // end namespace

#endif //KGL_RESOURCE_BASE_H