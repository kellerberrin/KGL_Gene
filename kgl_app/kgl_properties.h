//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H


#include "kgl_runtime.h"
#include "kel_property_tree.h"
#include "kgl_genome_types.h"
#include "kgl_properties_resource.h"

#include <memory>
#include <string>
#include <map>
#include <set>


namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// High level object extracts application specific properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeProperties {

public:

  RuntimeProperties() : property_tree_ptr_(std::make_shared<PropertyTree>()) {}
  ~RuntimeProperties() = default;

  [[nodiscard]] bool readProperties(const std::string& properties_file);

  void setWorkDirectory(const std::string& work_directory) { work_directory_ = work_directory; }

  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  [[nodiscard]] ActivePackageVector getActivePackages() const;

  [[nodiscard]] RuntimePackageMap getPackageMap() const;

  [[nodiscard]]  RuntimeAnalysisMap getAnalysisMap() const;

  [[nodiscard]]  ResourceDefinitions getRuntimeResources() const {

    ResourceProperties resource_properties(property_tree_ptr_, work_directory_);
    return resource_properties.getRuntimeResources();

  }

  [[nodiscard]] RuntimeDataFileMap getDataFiles() const;

  [[nodiscard]] ContigAliasMap getContigAlias() const;

  [[nodiscard]] VariantEvidenceMap getEvidenceMap() const;

  [[nodiscard]] ActiveParameterList getParameterMap() const;

private:

  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'
  std::shared_ptr<PropertyTree> property_tree_ptr_;   // The aggregated and parsed XML property tree.

  // Node categories.
  constexpr static const char DOT_[] = ".";
  constexpr static const char RUNTIME_ROOT_[] = "runTime";
  constexpr static const char ACTIVE_[] = "active";
  constexpr static const char VALUE_[] = "value";

  // Active Package Runtime categories.
  constexpr static const char EXECUTE_LIST_[] = "executeList";
  // Package Runtime categories.
  constexpr static const char PACKAGE_LIST_[] = "packageList";
  constexpr static const char PACKAGE_[] = "package";
  constexpr static const char PACKAGE_IDENT_[] = "packageIdent";
  constexpr static const char PACKAGE_ANALYSIS_LIST_[] = "analysisList";
  constexpr static const char PACKAGE_RESOURCE_LIST_[] = "resourceList";
  constexpr static const char PACKAGE_ITERATION_[] = "iteration";
  constexpr static const char PACKAGE_ITERATION_LIST_[] = "iterationList";
  // Analysis Runtime categories.
  constexpr static const char ANALYSIS_LIST_[] = "analysisList";
  constexpr static const char ANALYSIS_[] = "analysis";
  constexpr static const char ANALYSIS_IDENT_[] = "analysisIdent";
  // Parameter Runtime categories.
  constexpr static const char PARAMETER_LIST_[] = "parameterList";
  constexpr static const char PARAMETER_BLOCK_[] = "parameterBlock";
  constexpr static const char PARAMETER_NAME_[] = "parameterName";
  constexpr static const char PARAMETER_VECTOR_[] = "parameterVector";
  constexpr static const char PARAMETER_[] = "parameter";
  constexpr static const char PARAMETER_IDENT_[] = "parameterIdent";
  constexpr static const char PARAMETER_VALUE_[] = "parameterValue";
  constexpr static const char PARAMETER_RUNTIME_[] = "parameterRuntime";
  // Data File Runtime categories.
  constexpr static const char DATA_FILE_LIST_[] = "dataFileList";
  constexpr static const char DATA_FILE_IDENT_[] = "dataFileIdent";
  constexpr static const char DATA_FILE_NAME_[] = "dataFileName";
  constexpr static const char DATA_PARSER_TYPE_[] = "dataFileType";
  // general purpose data file without type specific fields.
  constexpr static const char GENERAL_DATA_FILE_TYPE_[] = "generalFile";
  // VCF data file specific fields.
  constexpr static const char VCF_DATA_FILE_TYPE_[] = "vcfFile";
  constexpr static const char VCF_FILE_GENOME_[] =  "vcfGenome";
  constexpr static const char VCF_INFO_EVIDENCE_[] =  "vcfInfo";
  // VCF Info Evidence categories.
  constexpr static const char EVIDENCE_LIST_[] = "evidenceList";
  constexpr static const char EVIDENCE_IDENT_[] = "evidenceIdent";
  constexpr static const char EVIDENCE_INFO_LIST_[] = "vcfInfoList";
  constexpr static const char EVIDENCE_INFO_ITEM_[] = "vcfInfoItem";
   // Contig/Chromosome Alias categories.
  constexpr static const char ALIAS_LIST_[] = "aliasList";
  constexpr static const char ALIAS_IDENT_[] = "ident";
  constexpr static const char ALIAS_TYPE_[] = "chromosomeType";
  constexpr static const char ALIAS_ENTRY_[] = "alias";

};


}   // end namespace



#endif //KGL_PROPERTIES_H
