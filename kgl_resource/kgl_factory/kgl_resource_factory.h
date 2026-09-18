//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_RESOURCE_FACTORY_H
#define KGL_RESOURCE_FACTORY_H


#include "kgl_runtime_resource.h"
#include "kel_property_tree.h"

#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <map>


namespace kellerberrin::genome {   //  organization::project level namespace


/// A two-phase resource factory. One subclass exists per resource type.
/// Phase 1 (parseXML): XML sub-tree -> ResourceParameters (lightweight definition).
/// Phase 2 (load): ResourceParameters -> shared_ptr<const ResourceBase> (heavyweight instance).
class ResourceFactory {

public:

  ResourceFactory() = default;
  virtual ~ResourceFactory() = default;

  /// The resource type string; matches the XML element name under <resourceList>
  /// and the resource's ResourceBase::resourceType().
  [[nodiscard]] virtual std::string_view resourceType() const = 0;

  /// Phase 1: build the immutable definition from a parsed XML sub-tree.
  /// Returns std::nullopt on recoverable parse error (error is logged).
  [[nodiscard]] virtual std::optional<ResourceParameters> parseXML(const PropertyTree& sub_tree,
                                                                  std::string_view work_directory) const = 0;

  /// Phase 2: construct the heavyweight resource from its definition.
  /// Calls ExecEnv::log().critical() on unrecoverable construction errors.
  [[nodiscard]] virtual std::shared_ptr<const ResourceBase> load(const ResourceParameters& parameters) const = 0;

};


/// Locate the factory for a given resource type. Returns nullptr for unknown types.
/// Uses a function-local static map with std::less<> for transparent string_view lookup.
[[nodiscard]] const ResourceFactory* findResourceFactory(std::string_view resource_type);


}   // end namespace


#endif //KGL_RESOURCE_FACTORY_H