//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_RUNTIME_CONFIG_H
#define KGL_RUNTIME_CONFIG_H


#include "kgl_properties.h"



namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Aggregation object to hold the various runtime arguments.
// A const reference of this object is passed to the analysis software.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RuntimeConfiguration {

public:

  RuntimeConfiguration(const RuntimeProperties &runtime_options, std::string work_directory)
      : active_packages_(runtime_options.getActivePackages()),
        contig_alias_(runtime_options.getContigAlias()),
        data_file_map_(runtime_options.getDataFiles()),
        resource_map_(runtime_options.getRuntimeResources()),
        analysis_map_(runtime_options.getAnalysisMap()),
        package_map_(runtime_options.getPackageMap()),
        evidence_map_(runtime_options.getEvidenceMap()),
        defined_parameters_(runtime_options.getParameterMap()),
        work_directory_(std::move(work_directory)) { verifyPackages(); }

  ~RuntimeConfiguration() = default;

  [[nodiscard]] const ActivePackageVector &activePackages() const { return active_packages_; }

  [[nodiscard]] const ContigAliasMap &contigAlias() const { return contig_alias_; }

  [[nodiscard]] const RuntimeDataFileMap &dataFileMap() const { return data_file_map_; }

  [[nodiscard]] const RuntimeResourceMap &resourceMap() const { return resource_map_; }

  [[nodiscard]] const RuntimeAnalysisMap &analysisMap() const { return analysis_map_; }

  [[nodiscard]] const RuntimePackageMap &runtimePackageMap() const { return package_map_; }

  [[nodiscard]] const VariantEvidenceMap &evidenceMap() const { return evidence_map_; }

  [[nodiscard]] const ActiveParameterList &activeParameterList() const { return defined_parameters_; }

  [[nodiscard]] const std::string &workDirectory() const { return work_directory_; }

private:

  // The Runtime information loaded from the XML config files.
  const ActivePackageVector active_packages_;
  const ContigAliasMap contig_alias_;
  const RuntimeDataFileMap data_file_map_;
  const RuntimeResourceMap resource_map_;
  const RuntimeAnalysisMap analysis_map_;
  const RuntimePackageMap package_map_;
  const VariantEvidenceMap evidence_map_;
  const ActiveParameterList defined_parameters_;
  const std::string work_directory_;

  // Check the integrity of all the XML information.
  void verifyPackages() const;

};





} // namespace



#endif //KGL_RUNTIME_CONFIG_H
