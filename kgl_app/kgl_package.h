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

  ExecutePackage( const std::string& options_file,
                  const std::string& option_file_out,
                  const std::string& parsed_option_out ,
                  const std::string& reference_directory)
                  : runtime_properties_( options_file, option_file_out, parsed_option_out, reference_directory),
                    runtime_config_(runtime_properties_, reference_directory),
                    package_manager_(runtime_config_),
                    package_analysis_(runtime_config_) {}

  ~ExecutePackage() = default;

  /// Executes all the application logic.
  void executeActive() const;

private:

  const RuntimeProperties    runtime_properties_;
  const RuntimeConfiguration runtime_config_;
  const PackageManager       package_manager_;
  const PackageAnalysis      package_analysis_;

};


} // namespace


#endif //KGL_PACKAGE_H