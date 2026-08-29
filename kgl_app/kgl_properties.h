//
// Created by kellerberrin on 11/11/18.
//

#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H


#include "kgl_runtime.h"
#include "kgl_runtime_resource.h"
#include "kel_property_tree.h"
#include "kgl_genome_types.h"

#include <memory>
#include <optional>
#include <string>
#include <utility>


namespace kellerberrin::genome {   //  organization::project level namespace


/// High level object that reads and caches application-specific properties from "runtime.xml".
/// Each accessor delegates to a dedicated xml:: section parser.
class RuntimeProperties {

public:

  RuntimeProperties() : property_tree_ptr_(std::make_shared<PropertyTree>()) {}
  ~RuntimeProperties() = default;
  RuntimeProperties(const RuntimeProperties&) = delete;
  RuntimeProperties& operator=(const RuntimeProperties&) = delete;
  RuntimeProperties(RuntimeProperties&&) = delete;
  RuntimeProperties& operator=(RuntimeProperties&&) = delete;

  /// Reads and parses the properties from the specified XML file.
  [[nodiscard]] bool readProperties(const std::string& properties_file);

  /// Sets the work directory used to resolve relative file paths.
  void setWorkDirectory(const std::string& work_directory) { work_directory_ = work_directory; }
  /// Returns the work directory.
  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  /// Returns the list of active packages to be executed at runtime.
  [[nodiscard]] const ActivePackageVector& getActivePackages() const;
  /// Returns the map of package identifiers to RuntimePackage objects.
  [[nodiscard]] const RuntimePackageMap&   getPackageMap() const;
  /// Returns the map of analysis identifiers to RuntimeAnalysis objects.
  [[nodiscard]] const RuntimeAnalysisMap&  getAnalysisMap() const;
  /// Returns the resource definitions parsed from the runtime XML.
  [[nodiscard]] const ResourceDefinitions& getRuntimeResources() const;
  /// Returns the map of data file identifiers to file info objects.
  [[nodiscard]] const RuntimeDataFileMap&  getDataFiles() const;
  /// Returns the contig alias map for chromosome identifier lookup.
  [[nodiscard]] const ContigAliasMap&      getContigAlias() const;
  /// Returns the variant evidence map for VCF INFO field specifications.
  [[nodiscard]] const VariantEvidenceMap&  getEvidenceMap() const;
  /// Returns the active parameter list for analysis configuration.
  [[nodiscard]] const ActiveParameterList& getParameterMap() const;

private:

  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'
  std::shared_ptr<PropertyTree> property_tree_ptr_;   // The aggregated and parsed XML property tree.

  // Each section is parsed at most once because RuntimeProperties is logically read-only after readProperties().
  template <class T, class F>
  [[nodiscard]] const T& cache(std::optional<T>& cache_slot, F&& parse) const {

    if (not cache_slot) cache_slot = std::forward<F>(parse)();
    return *cache_slot;

  }

  mutable std::optional<ActivePackageVector> active_packages_cache_;
  mutable std::optional<RuntimePackageMap>   package_map_cache_;
  mutable std::optional<RuntimeAnalysisMap>  analysis_map_cache_;
  mutable std::optional<ResourceDefinitions> runtime_resources_cache_;
  mutable std::optional<RuntimeDataFileMap>  data_files_cache_;
  mutable std::optional<ContigAliasMap>      contig_alias_cache_;
  mutable std::optional<VariantEvidenceMap>  evidence_map_cache_;
  mutable std::optional<ActiveParameterList> parameter_map_cache_;

};


} // namespace


#endif //KGL_PROPERTIES_H