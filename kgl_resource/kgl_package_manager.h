//
// Created by kellerberrin on 26/08/26.
//

#ifndef KGL_PACKAGE_MANAGER_H
#define KGL_PACKAGE_MANAGER_H


#include "kel_exec_env.h"
#include "kgl_runtime_config.h"
#include "kgl_runtime_resource.h"
#include "kgl_variant_db_population.h"

#include <memory>
#include <string>


namespace kellerberrin::genome {   //  organization::project level namespace


/// Loads package resources and reads iterative data files.
/// Separates resource/data mechanics from execution orchestration.
class PackageManager {

public:

  explicit PackageManager(const RuntimeConfiguration& runtime_config) : runtime_config_(runtime_config) {}
  ~PackageManager() = default;

  /// Load all resources requested by a package using the ResourceFactory registry.
  [[nodiscard]] std::shared_ptr<const AnalysisResources> loadResources(const RuntimePackage& package) const;

  /// Read a single data file identified by its package iteration list.
  [[nodiscard]] std::shared_ptr<DataDB> readDataFile(const RuntimePackage& package,
                                                     const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                                     const std::string& data_file_ident) const;

private:

  const RuntimeConfiguration& runtime_config_;

};


}   // end namespace


#endif //KGL_PACKAGE_MANAGER_H