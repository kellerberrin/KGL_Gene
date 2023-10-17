//
// Created by kellerberrin on 21/6/20.
//

#ifndef KGL_ANALYSIS_VIRTUAL_H
#define KGL_ANALYSIS_VIRTUAL_H



#include "kgl_runtime.h"
#include "kgl_runtime_resource.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace



// Virtual interface for the analysis class
class VirtualAnalysis {

public:


  VirtualAnalysis() = default;
  virtual ~VirtualAnalysis() = default;

  // Functions redefined in super classes
  // The ident must match the ident used in the package XML.
  [[nodiscard]] virtual std::string ident() const = 0;

  // Setup the analytics to process VCF data.
  [[nodiscard]] virtual bool initializeAnalysis( const std::string& work_directory,
                                                 const ActiveParameterList& named_parameters,
                                                 const std::shared_ptr<const AnalysisResources>& resource_ptr) = 0;

  // Perform the genetic analysis per VCF file
  [[nodiscard]] virtual bool fileReadAnalysis(std::shared_ptr<const DataDB> data_object_ptr) = 0;

  // Perform the genetic analysis per iteration
  [[nodiscard]] virtual bool iterationAnalysis() = 0;

  // All VCF data has been presented, finalize analysis and write results
  [[nodiscard]] virtual bool finalizeAnalysis() = 0;

  // The key is the analysis ident the value is the factory function to create the corresponding analysis object.
  // Note the key value must match the corresponding analysis ident in the definition XML file.
  using AnalysisFactoryMap = std::map<std::string, std::function<std::unique_ptr<VirtualAnalysis>(void)>>;
  // The map to hold the ident, factory function pairs.
  // This is used by PackageAnalysis to dynamically create specified analysis objects.
  const static AnalysisFactoryMap analysis_factory_map_;

private:

};


}

#endif //KGL_KGL_ANALYSIS_VIRT_H
