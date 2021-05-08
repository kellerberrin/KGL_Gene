//
// Created by kellerberrin on 1/5/20.
//

#include "kgl_package.h"
#include "kgl_variant_factory_parsers.h"
#include "kgl_ontology_database.h"

namespace kgl = kellerberrin::genome;
namespace kol = kellerberrin::ontology;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::ExecutePackage::executeActive() const {

  // for all active packages.
  for (auto const& active_package : runtime_config_.activePackages()) {

    // Get package definition
    auto result = runtime_config_.runtimePackageMap().find(active_package.packageIdentifier());

    if (result == runtime_config_.runtimePackageMap().end()) {

      ExecEnv::log().error( "ExecutePackage::executeActive, Active package :{} could not find matching package definition",
                            active_package.packageIdentifier());

      continue; // execute next active package.

    }

    auto [package_ident, package] = *result;

    // Get reference genomes.
    ExecEnv::log().info("Load Runtime Resources for Package: {}", package_ident);
    std::shared_ptr<const AnalysisResources> resource_ptr = loadRuntimeResources(package);

    // Setup the analytics
    if (not package_analysis_.initializeAnalysis(package, resource_ptr)) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem initializing Analysis for Package: {}", package_ident);

    }

    // Iterate through the VCF files and update analytics.
    for (auto const& iterative_files : package.iterativeFileList()) {

      for (auto const& data_file : iterative_files) {

        std::shared_ptr<DataDB> data_ptr = readDataFile(package, resource_ptr, data_file);

        if (not package_analysis_.fileReadAnalysis(data_ptr)) {

          ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Read File Analysis for Package: {}", package_ident);

        }

      }

      if (not package_analysis_.iterationAnalysis()) {

        ExecEnv::log().error("ExecutePackage::executeActive, Problem performing Analysis for Package: {}", package_ident);

      }

    }

    // Complete and write the analytics.
    if (not package_analysis_.finalizeAnalysis()) {

      ExecEnv::log().error("ExecutePackage::executeActive, Problem finalizing Analysis for Package: {}", package_ident);

    }

  }

}


std::shared_ptr<kgl::DataDB>
kgl::ExecutePackage::readDataFile(const RuntimePackage& package,
                                  const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                  const std::string& data_file) const {


  ExecEnv::log().info("Package: {}, Data file ident: {}", package.packageIdentifier(), data_file);

  auto result = runtime_config_.dataFileMap().find(data_file);
  if (result == runtime_config_.dataFileMap().end()) {

    ExecEnv::log().critical("ExecutePackage::readDataFile, Package: {}, data file ident: {}, not defined",
                             package.packageIdentifier(), data_file);

  }

  auto [file_ident, file_info_ptr] = *result;

  // Selects the appropriate parser and returns a base class data object.
  std::shared_ptr<kgl::DataDB> data_ptr = ParserSelection::parseData(resource_ptr, file_info_ptr, runtime_config_.evidenceMap(), runtime_config_.contigAlias());

  return data_ptr;

}

