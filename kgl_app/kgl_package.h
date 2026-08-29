//
// Created by kellerberrin on 1/5/20.
//

#ifndef KGL_PACKAGE_H
#define KGL_PACKAGE_H

#include "kel_exec_env.h"
#include "kgl_runtime_config.h"
#include "kgl_package_analysis.h"
#include "kgl_package_manager.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/// Executes the analysis software and provides the software with data files and genome resources.
class ExecutePackage {

public:

  ExecutePackage(const RuntimeProperties& runtime_options, const std::string& work_directory)
                  : runtime_config_(runtime_options, work_directory),
                    package_manager_(runtime_config_),
                    package_analysis_(runtime_config_) {}

  ~ExecutePackage() = default;

  /// Executes all the application logic.
  void executeActive() const;

private:

  const RuntimeConfiguration runtime_config_;
  const PackageManager       package_manager_;
  const PackageAnalysis      package_analysis_;

};


} // namespace


#endif //KGL_PACKAGE_H