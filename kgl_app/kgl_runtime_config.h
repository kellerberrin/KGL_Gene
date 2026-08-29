//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_RUNTIME_CONFIG_H
#define KGL_RUNTIME_CONFIG_H


#include "kgl_properties.h"



namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Validated view over RuntimeProperties.
//
// Construction eagerly calls all RuntimeProperties accessors (triggering lazy parse)
// and runs verifyPackages() to validate all cross-references before any work begins.
// After construction, the configuration is a read-only view over the RuntimeProperties
// lazy caches — no data is copied.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Validated read-only view over RuntimeProperties.
/// Construction eagerly calls all RuntimeProperties accessors and runs verifyPackages() to validate cross-references.
class RuntimeConfiguration {

public:

  RuntimeConfiguration(const RuntimeProperties &runtime_options, std::string work_directory)
      : runtime_options_(runtime_options),
        work_directory_(std::move(work_directory)) { verifyPackages(); }

  ~RuntimeConfiguration() = default;

  /// Returns the list of active packages to be executed.
  [[nodiscard]] const ActivePackageVector &activePackages() const { return runtime_options_.getActivePackages(); }

  /// Returns the contig alias map for chromosome identifier lookup.
  [[nodiscard]] const ContigAliasMap &contigAlias() const { return runtime_options_.getContigAlias(); }

  /// Returns the map of data file identifiers to file info objects.
  [[nodiscard]] const RuntimeDataFileMap &dataFileMap() const { return runtime_options_.getDataFiles(); }

  /// Returns the resource definitions parsed from the runtime XML.
  [[nodiscard]] const ResourceDefinitions &resourceDefMap() const { return runtime_options_.getRuntimeResources(); }

  /// Returns the map of analysis identifiers to RuntimeAnalysis objects.
  [[nodiscard]] const RuntimeAnalysisMap &analysisMap() const { return runtime_options_.getAnalysisMap(); }

  /// Returns the map of package identifiers to RuntimePackage objects.
  [[nodiscard]] const RuntimePackageMap &runtimePackageMap() const { return runtime_options_.getPackageMap(); }

  /// Returns the variant evidence map for VCF INFO field specifications.
  [[nodiscard]] const VariantEvidenceMap &evidenceMap() const { return runtime_options_.getEvidenceMap(); }

  /// Returns the active parameter list for analysis configuration.
  [[nodiscard]] const ActiveParameterList &activeParameterList() const { return runtime_options_.getParameterMap(); }

  /// Returns the work directory.
  [[nodiscard]] const std::string &workDirectory() const { return work_directory_; }

private:

  const RuntimeProperties& runtime_options_;
  const std::string work_directory_;

  // Check the integrity of all the XML information.
  void verifyPackages() const;

};





} // namespace



#endif //KGL_RUNTIME_CONFIG_H