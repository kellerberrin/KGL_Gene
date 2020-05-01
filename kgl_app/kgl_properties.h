//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H


#include "kgl_runtime.h"
#include "kel_property_tree.h"
#include "kgl_genome_types.h"

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

  RuntimeProperties() = default;
  ~RuntimeProperties() = default;

  [[nodiscard]] bool readProperties(const std::string& properties_file);

  void setWorkDirectory(const std::string& work_directory) { work_directory_ = work_directory; }

  [[nodiscard]] const std::string& workDirectory() const { return work_directory_; }

  [[nodiscard]] RuntimePackageMap getPackageMap() const;

  [[nodiscard]]  RuntimeAnalysisMap getAnalysisMap() const;

  [[nodiscard]]  RuntimeGenomeDatabaseMap getGenomeReferenceMap() const;

  [[nodiscard]] RuntimeVCFFileMap getVCFFiles() const;

  [[nodiscard]] ContigAliasMap getContigAlias() const;

////////////////////////////////////////////////////////////////////////////////////////////////////
// Legacy Code.

  [[nodiscard]] bool getMixtureFile(std::string& mixture_file) const;

  [[nodiscard]] bool getPropertiesAuxFile(std::string &aux_file) const;


  void getGenomeDBFiles(const std::string& organism,
                        std::string& fasta_file,
                        std::string& gff_file,
                        std::string& gaf_file,
                        std::string& tranlation_table) const;

  [[nodiscard]] bool getGenomeAuxFiles(const std::string& organism, std::vector<AuxFileInfo>& auxfiles) const;


  [[nodiscard]] bool getActiveGenomes(std::vector<std::string>& genome_list) const;

private:

  // Node categories.
  constexpr static const char* DOT_ = ".";
  constexpr static const char* RUNTIME_ROOT_ = "runTime";
  constexpr static const char* HELP_ = "help";
  constexpr static const char* ACTIVE_ = "active";
  constexpr static const char* VALUE_ = "value";
  // Package Runtime categories.
  constexpr static const char* PACKAGE_LIST_ = "packageList";
  constexpr static const char* PACKAGE_ = "package";
  constexpr static const char* PACKAGE_IDENT_ = "packageIdent";
  constexpr static const char* PACKAGE_ANALYSIS_LIST_ = "analysisList";
  constexpr static const char* PACKAGE_GENOME_LIST_ = "genomeList";
  constexpr static const char* PACKAGE_LOAD_LIST_ = "loadList";
  constexpr static const char* PACKAGE_ITERATION_LIST_ = "iterationList";
  // Analysis Runtime categories.
  constexpr static const char* ANALYSIS_LIST_ = "analysisList";
  constexpr static const char* ANALYSIS_ = "analysis";
  constexpr static const char* ANALYSIS_IDENT_ = "analysisIdent";
  constexpr static const char* PARAMETER_LIST_ = "analysisIdent";
  constexpr static const char* PARAMETER_ = "parameter";
  constexpr static const char* PARAMETER_IDENT_ = "parameterIdent";
  constexpr static const char* PARAMETER_VALUE_ = "parameterIdent";
  // VCF File Runtime categories.
  constexpr static const char* VCF_LIST_ = "vcfList";
  constexpr static const char* VCF_IDENT_ = "vcfIdent";
  constexpr static const char* VCF_FILE_ = "vcfFile";
  constexpr static const char* VCF_FILE_NAME_ = "vcfFileName";
  constexpr static const char* VCF_PARSER_TYPE_ = "vcfParser";
  constexpr static const char* VCF_FILE_GENOME_ =  "vcfGenome";
  constexpr static const char* VCF_FILE_PLOIDY_ =  "vcfPloidy";
  // Genome Database categories.
  constexpr static const char* GENOME_LIST_ = "genomeList";
  constexpr static const char* GENOME_DATABASE_ = "genomeDatabase";
  constexpr static const char* GENOME_IDENT_ = "genomeIdent";
  constexpr static const char* FASTA_FILE_ = "fastaFile";
  constexpr static const char* GFF_FILE_ = "gffFile";
  constexpr static const char* GAF_FILE_ = "gafFile";   // optional
  constexpr static const char* TRANSLATION_TABLE_ = "translationTable";
  // Contig/Chromosome Alias categories.
  constexpr static const char* ALIAS_LIST_ = "contigAlias";
  constexpr static const char* ALIAS_IDENT_ = "ident";
  constexpr static const char* ALIAS_ENTRY_ = "alias";
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Legacy
  constexpr static const char* FILE_LIST_ = "fileList";
  constexpr static const char* VCF_PLOIDY_ = "vcfPloidy";
  constexpr static const char* MIXTURE_FILE_ = "mixtureFile";
  constexpr static const char* AUX_FILE_ = "auxFile";
  constexpr static const char* AUX_GENOME_FILE_ = "auxGenomeFile";
  // Defaults
  constexpr static const size_t DEFAULT_PLOIDY_ = 2;

  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'
  PropertyTree property_tree_;   // All the option XML files.

};


}   // end namespace



#endif //KGL_PROPERTIES_H
