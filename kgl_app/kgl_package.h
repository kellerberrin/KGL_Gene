//
// Created by kellerberrin on 1/5/20.
//

#ifndef KGL_PACKAGE_H
#define KGL_PACKAGE_H

#include "kel_exec_env.h"
#include "kgl_runtime_config.h"
#include "kgl_runtime_resource.h"
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

  // Load the package resources.
  [[nodiscard]] std::shared_ptr<const AnalysisResources> loadRuntimeResources(const RuntimePackage& package) const;
  void loadGenomeResource(const std::string& resource_type,  const std::string& genome_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadOntologyResource(const std::string& resource_type, const std::string& ontology_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadHsGeneNomenclatureResource(const std::string& resource_type, const std::string& nomenclature_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadHsGenomeGenealogyResource(const std::string& resource_type, const std::string& genealogy_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadHsGenomeAuxResource(const std::string& resource_type, const std::string& genome_aux_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadCitationResource(const std::string& resource_type, const std::string& citation_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadEntrezGeneResource(const std::string& resource_type, const std::string& entrez_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadPubmedAPIResource(const std::string& resource_type, const std::string& api_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadPf7SampleResource(const std::string& resource_type, const std::string& Pf7_sample_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadPf7FwsResource(const std::string& resource_type, const std::string& Pf7_Fws_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;
  void loadPf7DistanceResource(const std::string& resource_type, const std::string& Pf7_Fws_ident, const std::shared_ptr<AnalysisResources>& resource_ptr) const;

  // Load a specified data file and return a base pointer (DataDB) to the file.
  [[nodiscard]] std::shared_ptr<DataDB> readDataFile(const RuntimePackage& package,
                                                     const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                                     const std::string& data_file) const;



};




} // namespace

#endif //KGL_PACKAGE_H
