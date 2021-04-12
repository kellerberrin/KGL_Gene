//
// Created by kellerberrin on 1/5/20.
//

#ifndef KGL_PACKAGE_H
#define KGL_PACKAGE_H

#include "kel_exec_env.h"
#include "kgl_runtime_config.h"
#include "kgl_genome_collection.h"
#include "kgl_package_analysis.h"

namespace kellerberrin::genome {   //  organization::project level namespace



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Executes the analysis software and provides the software with data files and genome resources.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ExecutePackage {

public:

  ExecutePackage( const RuntimeProperties& runtime_options, const std::string& work_directory)
                  : runtime_config_(runtime_options, work_directory),
                    package_analysis_(runtime_config_) {}

  ~ExecutePackage() = default;

  // Executes all the application logic.
  void executeActive() const;

private:

  // The Runtime information loaded from the XML config files.
  const RuntimeConfiguration runtime_config_;
  // The analysis management object.
  const PackageAnalysis package_analysis_;

  // Load the reference genomes.
  [[nodiscard]] std::unique_ptr<GenomeCollection> loadRuntimeResources(const RuntimePackage& package) const;

  // Load a specified data file and return a base pointer (DataDB) to the file.
  [[nodiscard]] std::shared_ptr<DataDB> readDataFile(const RuntimePackage& package,
                                                     std::shared_ptr<const GenomeCollection> reference_genomes,
                                                     const std::string& data_file) const;



};




} // namespace

#endif //KGL_PACKAGE_H
