//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_runtime_config.h"



namespace kgl = kellerberrin::genome;



/// Verifies that all packages, analyses, resources, and data files referenced in the runtime configuration are defined.
void kgl::RuntimeConfiguration::verifyPackages() const {

  // List the active packages.
  for (auto const& active : runtime_options_.getActivePackages()) {

    ExecEnv::log().info("Sequentially Executing Active Package: {}", active.packageIdentifier());

  }

  // for all packages.
  for (auto const& [package_ident, package] : runtime_options_.getPackageMap()) {

    // Confirm that requested analytics are defined.
    for (auto const& analysis_ident : package.analysisList()) {

      if (not runtime_options_.getAnalysisMap().contains(analysis_ident)) {

        ExecEnv::log().error("RuntimeConfiguration::verifyPackages, Package: {}, Analysis: {}, not defined", package_ident, analysis_ident);
        for (auto const& [ident, runtime] : runtime_options_.getAnalysisMap()) {

          ExecEnv::log().warn("RuntimeConfiguration::verifyPackages, Available Analysis: {}", ident);

        }

        ExecEnv::log().critical("RuntimeConfiguration::verifyPackages, Package: {}, Analysis: {}, not defined, unrecoverable.", package_ident, analysis_ident);

      }

      // If the analysis exists, check that any active named parameter blocks also exist.
      auto const& [analysis_id, analysis_obj] = *(runtime_options_.getAnalysisMap().find(analysis_ident));
      for (auto const& param_name : analysis_obj.parameterMap()) {

        if (not runtime_options_.getParameterMap().getMap().contains(param_name)) {

          ExecEnv::log().critical("RuntimeConfiguration::verifyPackages, Package: {}, Analysis: {}, Named Parameter Block: {} not defined",
                                  package_ident, analysis_ident, param_name);

        }

      }

    }

    //Confirm that all resources exist.
    for (auto const& [resource_type, resource_ident] : package.resourceList()) {

      if (not runtime_options_.getRuntimeResources().retrieve(resource_type, resource_ident)) {

        for (auto const& [ident, resource_parameters] : runtime_options_.getRuntimeResources().getMap()) {

          ExecEnv::log().info("RuntimeConfiguration::verifyPackages, Package: {}, Resource Type: {}, Ident: {}",
                              package_ident, resource_parameters.resourceType(), resource_parameters.resourceIdent());

        }

        ExecEnv::log().critical("RuntimeConfiguration::verifyPackages, Package: {}, Runtime Resource Type: '{}', Ident: '{}' not defined",
                                package_ident, resource_type, resource_ident);

      }

    }

    //confirm that all iterative load files exist.
    // Note that iterativeFileList() returns a nested vector, std::vector<std::vector<std::string>>
    for (auto const& vcf_file_vector : package.iterativeFileList()) {

      for (auto const& vcf_file_ident : vcf_file_vector) {

        if (not runtime_options_.getDataFiles().contains(vcf_file_ident)) {

          ExecEnv::log().critical("RuntimeConfiguration::verifyPackages, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

        }

      }

    }

    ExecEnv::log().info("Package: {}, All Reference Genomes, data files and analysis types are defined.", package_ident);

  }

}