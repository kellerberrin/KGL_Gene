//
// Created by kellerberrin on 15/1/21.
//

#include "kgl_runtime_config.h"



namespace kgl = kellerberrin::genome;



void kgl::RuntimeConfiguration::verifyPackages() const {

  // List the active packages.
  for (auto const& active : active_packages_) {

    ExecEnv::log().info("Sequentially Executing Active Package: {}", active.packageIdentifier());

  }

  // for all packages.
  for (auto const& [package_ident, package] : package_map_) {

    // Confirm that requested analytics are defined.
    for (auto const& analysis_ident : package.analysisList()) {

      if (not analysis_map_.contains(analysis_ident)) {

        ExecEnv::log().error("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, not defined", package_ident, analysis_ident);
        for (auto const& [ident, runtime] : analysis_map_) {

          ExecEnv::log().warn("ExecutePackage::verifyPackage, Available Analysis: {}", ident);

        }

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, not defined, unrecoverable.", package_ident, analysis_ident);

      }

      // If the analysis exists, check that any active named parameter blocks also exist.
      auto const& [analysis_id, analysis_obj] = *(analysis_map_.find(analysis_ident));
      for (auto const& param_name : analysis_obj.parameterMap()) {

        if (not defined_parameters_.getMap().contains(param_name)) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, Named Parameter Block: {} not defined",
                                  package_ident, analysis_ident, param_name);

        }

      }

    }

    //confirm that all reference genomes exist
    for (auto const& [resource_type, resource_ident] : package.resourceDatabaseList()) {

      if (not resource_map_.contains(resource_ident)) {

        for (auto const& [id, type] : resource_map_) {

          ExecEnv::log().info("ExecutePackage::verifyPackage, Package: {}, Resource map content: {}", package_ident, id);

        }

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Runtime Resource: {}, not defined", package_ident, resource_ident);

      }

    }

    //confirm that all iterative load files exist.
    // Note that iterativeFileList() returns a nested vector, std::vector<std::vector<std::string>>
    for (auto const& vcf_file_vector : package.iterativeFileList()) {

      for (auto const& vcf_file_ident : vcf_file_vector) {

        if (not data_file_map_.contains(vcf_file_ident)) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

        }

      }

    }

    ExecEnv::log().info("Package: {}, All Reference Genomes, data files and analysis types are defined.", package_ident);

  }

}

