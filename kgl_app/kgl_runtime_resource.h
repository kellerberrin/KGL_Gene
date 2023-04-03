//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_RESOURCE_DB_H
#define KGL_RESOURCE_DB_H

#include "kel_exec_env.h"

#include <optional>
#include <memory>
#include <string>
#include <vector>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ResourceBase inherited by all resources.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ResourceBase {

public:

  ResourceBase(std::string resource_type, std::string resource_ident)
      : resource_type_(std::move(resource_type)), resource_ident_(std::move(resource_ident)) {}
  virtual ~ResourceBase() = default;

  [[nodiscard]] const std::string& resourceType() const { return resource_type_; }
  [[nodiscard]] const std::string& resourceIdent() const { return resource_ident_; }

private:

  const std::string resource_type_;
  const std::string resource_ident_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simple container object to hold resources as they are passed to analysis packages.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using ResourceMap = std::multimap<std::string, std::shared_ptr<const ResourceBase>>;
class AnalysisResources {

public:

  AnalysisResources() = default;
  ~AnalysisResources() = default;

  void addResource(const std::shared_ptr<const ResourceBase>& resource_ptr) {

    resource_map_.emplace(resource_ptr->resourceType(), resource_ptr);

  }



  [[nodiscard]] std::vector<std::shared_ptr<const ResourceBase>> getResources( const std::string& resource_type,
                                                                               const std::string& resource_ident = "") const;

  [[nodiscard]] const ResourceMap& getMap() const { return resource_map_; }


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


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects encode XML resource definitions and pass then to the resource constructor.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ResourceParameterMap = std::map<std::string, std::string>;
class ResourceParameters {

public:

  ResourceParameters(std::string resource_type, std::string resource_id) : resource_type_(std::move(resource_type)),
                                                                           resource_id_(std::move(resource_id)) {}
  ~ResourceParameters() = default;

  [[nodiscard]] std::optional<const std::string> getParameter(const std::string& parameter_key) const;
  void setParameter(const std::string& parameter_key, const std::string& parameter_value);

  [[nodiscard]] const std::string& resourceType() const { return resource_type_; }
  [[nodiscard]] const std::string& resourceIdent() const { return resource_id_; }

private:

  std::string resource_type_;
  std::string resource_id_;
  ResourceParameterMap parameter_map_;

};



using ResourceDefinitionsMap = std::multimap<std::string, ResourceParameters>;
class ResourceDefinitions {

public:

  ResourceDefinitions() = default;
  ResourceDefinitions(const ResourceDefinitions&) = default;
  ~ResourceDefinitions() = default;

  void insert(const std::pair<std::string, ResourceParameters>& resource_def) { resource_definition_map_.insert(resource_def); }
  [[nodiscard]] std::optional<ResourceParameters> retrieve(const std::string& resource_type, const std::string& resource_ident) const;

  [[nodiscard]] const ResourceDefinitionsMap& getMap() const { return resource_definition_map_; }

private:

  ResourceDefinitionsMap resource_definition_map_;

};



}   // end namespace

#endif //KGL_RESOURCE_DB_H
