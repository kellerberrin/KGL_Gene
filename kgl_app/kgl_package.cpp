//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"


namespace kgl = kellerberrin::genome;


/// Sequentially executes all active packages defined in the runtime configuration.
void kgl::ExecutePackage::executeActive() const {

  for (auto const& active_package : runtime_config_.activePackages()) {

    auto result = runtime_config_.runtimePackageMap().find(active_package.packageIdentifier());
    if (result == runtime_config_.runtimePackageMap().end()) {

      ExecEnv::log().error("ExecutePackage::executeActive, Active package :{} could not find matching package definition",
                           active_package.packageIdentifier());
      continue;

    }

    auto [package_ident, package] = *result;

    ExecEnv::log().info("Load Runtime Resources for Package: {}", package_ident);
    std::shared_ptr<const AnalysisResources> resource_ptr = package_manager_.loadResources(package);

    if (not package_analysis_.initializeAnalysis(package, resource_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem initializing Analysis for Package: {}", package_ident);

    }

    for (auto const& iterative_files : package.iterativeFileList()) {

      for (auto const& data_file : iterative_files) {

        std::shared_ptr<DataDB> data_ptr = package_manager_.readDataFile(package, resource_ptr, data_file);
        if (not package_analysis_.fileReadAnalysis(data_ptr)) {

          ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Read File Analysis for Package: {}", package_ident);

        }

      }

      if (not package_analysis_.iterationAnalysis()) {

        ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Analysis for Package: {}", package_ident);

      }

    }

    if (not package_analysis_.finalizeAnalysis()) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem finalizing Analysis for Package: {}", package_ident);

    }

  }

}