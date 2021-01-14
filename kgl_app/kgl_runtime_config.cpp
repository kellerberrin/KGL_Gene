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

      auto result = analysis_map_.find(analysis_ident);
      if (result == analysis_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, not defined", package_ident, analysis_ident);

      }

      // If the analysis exists, check that any active named parameter blocks also exist.
      auto const& [analysis_id, analysis_obj] = *result;

      for (auto const& param_name : analysis_obj.parameterMap()) {

        auto param_result = defined_parameters_.getMap().find(param_name);

        if (param_result ==  defined_parameters_.getMap().end()) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Analysis: {}, Named Parameter Block: {} not defined",
                                  package_ident, analysis_ident, param_name);

        }

      }

    }

    //confirm that all reference genomes exist
    for (auto const& genome_ident : package.genomeDatabaseList()) {

      auto result = genome_map_.find(genome_ident);
      if (result == genome_map_.end()) {

        ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Reference Genome: {}, not defined", package_ident, genome_ident);

      }

    }

    //confirm that all iterative load files exist.
    // Note that iterativeFileList() returns a nested vector, std::vector<std::vector<std::string>>
    for (auto const& vcf_file_vector : package.iterativeFileList()) {

      for (auto const& vcf_file_ident : vcf_file_vector) {

        auto result = data_file_map_.find(vcf_file_ident);
        if (result == data_file_map_.end()) {

          ExecEnv::log().critical("ExecutePackage::verifyPackage, Package: {}, Iterative load file: {}, not defined", package_ident, vcf_file_ident);

        }

      }

    }

    ExecEnv::log().info("Package: {}, All Reference Genomes, data files and analysis types are defined.", package_ident);

  }

}

