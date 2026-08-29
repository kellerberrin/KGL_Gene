//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_RESOURCE_PARAMS_H
#define KGL_RESOURCE_PARAMS_H

#include <optional>
#include <string>
#include <map>

namespace kellerberrin::genome {   //  organization level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// These objects encode XML resource definitions and pass them to the resource constructor.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ResourceParameterMap = std::map<std::string, std::string>;

/// Holds the parameters for a single resource definition parsed from XML.
class ResourceParameters {

public:

  ResourceParameters(std::string resource_type, std::string resource_id) : resource_type_(std::move(resource_type)),
                                                                           resource_id_(std::move(resource_id)) {}
  ~ResourceParameters() = default;

  /// Looks up a parameter value by key. Returns std::nullopt if not found.
  [[nodiscard]] std::optional<const std::string> getParameter(const std::string& parameter_key) const;
  /// Sets a parameter value. Warns on duplicate key.
  void setParameter(const std::string& parameter_key, const std::string& parameter_value);

  /// Returns the resource type string.
  [[nodiscard]] const std::string& resourceType() const { return resource_type_; }
  /// Returns the resource identifier string.
  [[nodiscard]] const std::string& resourceIdent() const { return resource_id_; }

private:

  std::string resource_type_;
  std::string resource_id_;
  ResourceParameterMap parameter_map_;

};



using ResourceDefinitionsMap = std::multimap<std::string, ResourceParameters>;

/// A multimap of resource definitions keyed by resource type. Allows multiple resources of the same type with different identifiers.
class ResourceDefinitions {

public:

  ResourceDefinitions() = default;
  ResourceDefinitions(const ResourceDefinitions&) = default;
  ~ResourceDefinitions() = default;

  /// Inserts a resource definition into the map.
  void insert(const std::pair<std::string, ResourceParameters>& resource_def) { resource_definition_map_.insert(resource_def); }
  /// Looks up a resource definition by type and identifier. Returns std::nullopt if not found.
  [[nodiscard]] std::optional<ResourceParameters> retrieve(const std::string& resource_type, const std::string& resource_ident) const;

  /// Returns the underlying multimap of all resource definitions.
  [[nodiscard]] const ResourceDefinitionsMap& getMap() const { return resource_definition_map_; }

private:

  ResourceDefinitionsMap resource_definition_map_;

};


}   // end namespace

#endif //KGL_RESOURCE_PARAMS_H