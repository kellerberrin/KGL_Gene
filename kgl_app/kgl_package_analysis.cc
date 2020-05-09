//
// Created by kellerberrin on 4/5/20.
//

#include "kgl_package_analysis.h"


namespace kgl = kellerberrin::genome;


bool kgl::PackageAnalysis::initializeAnalysis( const RuntimePackage& package,
                                              std::shared_ptr<const GenomeCollection> reference_genomes) const {


  active_analysis_.clear();
  for (auto const& analysis_id : package.analysisList()) {

    bool found = false;
    for (auto const& registered_analysis : registered_analysis_) {

      if (not registered_analysis) {

        ExecEnv::log().error("Registered Analysis is a NULL pointer");
        continue;

      }

      if (registered_analysis->ident() == analysis_id) {

        found = true;

        std::unique_ptr<NullAnalysis> active_analysis = registered_analysis->factory();

        auto result = analysis_map_.find(analysis_id);

        // Analysis parameters found.
        if (result != analysis_map_.end()) {

          if (active_analysis->initializeAnalysis(work_directory_, result->second.parameterMap(), reference_genomes)) {

            active_analysis_.emplace_back(std::move(active_analysis), true);

          } else {

            active_analysis_.emplace_back(std::move(active_analysis), false);  // Register but disable.
            ExecEnv::log().error("Failed to Initialize Analysis: {}, further analysis discarded", analysis_id);

          }

        } else {

          ExecEnv::log().error("Could not find Runtime Parameters for Analysis: {}, further analysis discarded", analysis_id);

        }

        break;

      }

    }

    if (not found) {

      ExecEnv::log().error("Could not find Analysis: {}. Analysis must be registered in PackageAnalysis::PackageAnalysis", analysis_id);

    }

  }

  return true;
}


bool kgl::PackageAnalysis::iterateAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes,
                                           std::shared_ptr<const UnphasedPopulation> vcf_iterative_data) const {

  for (auto& [analysis, active] : active_analysis_) {

    ExecEnv::log().info("PackageAnalysis::iterateAnalysis; Updating Analysis: {}", analysis->ident());

    if (active) {

      if (not analysis->iterateAnalysis(reference_genomes, vcf_iterative_data)) {

        ExecEnv::log().error("PackageAnalysis::iterateAnalysis; Error Iteratively Updating Analysis: {}, disabled from further updates.", analysis->ident());
        active = false;
        return false;

      }

    } else {

      ExecEnv::log().warn("PackageAnalysis::iterateAnalysis; Analysis: {} registered an error code and was disabled.", analysis->ident());

    }

  }

  return true;

}


bool kgl::PackageAnalysis::finalizeAnalysis(std::shared_ptr<const GenomeCollection> reference_genomes) const {


  for (auto& [analysis, active] : active_analysis_) {

    if (active) {

      if (not analysis->finalizeAnalysis(reference_genomes)) {

        ExecEnv::log().error("Error Finalizing Updating Analysis: {}", analysis->ident());
        active = false;
        return false;

      }

    } else {

      ExecEnv::log().warn("PackageAnalysis::finalizeAnalysis; Analysis: {} registered an error code and was disabled.", analysis->ident());

    }

  }

  return true;

}
