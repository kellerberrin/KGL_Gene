//
// Created by kellerberrin on 4/5/20.
//

#include "kgl_package_analysis.h"


namespace kgl = kellerberrin::genome;


bool kgl::PackageAnalysis::initializeAnalysis( const RuntimePackage& package,
                                               const std::shared_ptr<const AnalysisResources>& resource_ptr) const {


  active_analysis_.clear();
  // There can be multiple analysis defined for a package.
  for (auto const& analysis_id : package.analysisList()) {

    // Find the corresponding factory function.
    auto find_iter = VirtualAnalysis::analysis_factory_map_.find(analysis_id);
    if (find_iter != VirtualAnalysis::analysis_factory_map_.end()) {

      auto [key_id, factory] = *find_iter;
      std::unique_ptr<VirtualAnalysis> analysis_ptr = factory();

      auto result = runtime_contig_.analysisMap().find(analysis_id);
      // Analysis parameters found.
      if (result != runtime_contig_.analysisMap().end()) {

        // Get the named parameter blocks for this analysis (if any).
        auto [runtime_id, analysis_runtime] = *result;
        auto defined_parameters = runtime_contig_.activeParameterList().createParameterList(analysis_runtime.parameterMap());

        // Initialize the analysis.
        if (analysis_ptr->initializeAnalysis(runtime_contig_.workDirectory(), defined_parameters, resource_ptr)) {

          active_analysis_.emplace_back(std::move(analysis_ptr), true);

        } else {

          active_analysis_.emplace_back(std::move(analysis_ptr), false);  // Register but disable.
          ExecEnv::log().error("PackageAnalysis::initializeAnalysis, Failed to Initialize Analysis: {}, further analysis discarded", analysis_id);

        }

      } else {

        ExecEnv::log().error("PackageAnalysis::initializeAnalysis, Could not find Runtime Parameters for Analysis: {}, further analysis discarded", analysis_id);

      }

    } else {

      ExecEnv::log().error("PackageAnalysis::initializeAnalysis, Could not find Analysis: {}. Analysis must be registered in PackageAnalysis::PackageAnalysis", analysis_id);

    }

  }

  return true;
}


bool kgl::PackageAnalysis::fileReadAnalysis(std::shared_ptr<const DataDB> file_data) const {

  for (auto& [analysis, active] : active_analysis_) {

    ExecEnv::log().info("File Read, Updating Analysis: {}", analysis->ident());

    if (active) {

      if (not analysis->fileReadAnalysis(file_data)) {

        ExecEnv::log().error("PackageAnalysis::fileReadAnalysis; Error Iteratively Updating Analysis: {}, disabled from further updates.", analysis->ident());
        active = false;
        return false;

      }

    } else {

      ExecEnv::log().warn("PackageAnalysis::fileReadAnalysis; Analysis: {} registered an error code and was disabled.", analysis->ident());

    }

  }

  ExecEnv::log().info("PackageAnalysis::fileReadAnalysis; Completed all file analysis");

  return true;

}


bool kgl::PackageAnalysis::iterationAnalysis() const {

  for (auto& [analysis, active] : active_analysis_) {

    ExecEnv::log().info("Updating Analysis: {}", analysis->ident());

    if (active) {

      if (not analysis->iterationAnalysis()) {

        ExecEnv::log().error("PackageAnalysis::iterationAnalysis; Error Iteratively Updating Analysis: {}, disabled from further updates.", analysis->ident());
        active = false;
        return false;

      }

    } else {

      ExecEnv::log().warn("PackageAnalysis::iterationAnalysis; Analysis: {} registered an error code and was disabled.", analysis->ident());

    }

  }

  return true;

}


bool kgl::PackageAnalysis::finalizeAnalysis() const {


  for (auto& [analysis, active] : active_analysis_) {

    if (active) {

      if (not analysis->finalizeAnalysis()) {

        ExecEnv::log().error("PackageAnalysis::initializeAnalysis, Error Finalizing Updating Analysis: {}", analysis->ident());
        active = false;
        return false;

      }

    } else {

      ExecEnv::log().warn("PackageAnalysis::finalizeAnalysis; Analysis: {} registered an error code and was disabled.", analysis->ident());

    }

  }

  return true;

}
