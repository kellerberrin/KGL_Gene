//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H

#include "kgl_exec_env.h"

#include <memory>
#include <string>
#include <set>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



// This object is a facade over the boost:: property tree object
// Properties may stored as XML or JSON. The actual format is determined at runtime
// by a file extension of "xml" or "json". Any other file extension is assumed to be
// formatted as XML. The choice of which of the two formats to use is at the
// discretion of the program user and is transparent to this object.
// The boost:: functionality is hidden using the PIMPL idiom.

class PropertyTree; // fwd
using SubPropertyTree = std::pair<std::string, PropertyTree>;

class PropertyTree {

public:

  class PropertyImpl;       // Forward declaration of the boost:: properties implementation class

  PropertyTree(); // defined in implementation file
  PropertyTree(const PropertyTree&);
  PropertyTree(const PropertyImpl&);

  ~PropertyTree(); // defined in implementation file

  bool readProperties(const std::string& properties_file);

  bool getProperty(const std::string& property_name, std::string& property) const;

  bool getPropertyVector(const std::string& property_name, std::vector<std::string>& property_vector) const;

  bool getPropertyTreeVector(const std::string& property_name, std::vector<SubPropertyTree>& property_tree_vector) const;

  bool getProperty(const std::string& property_name, size_t& property) const;

  bool getFileProperty(const std::string& property_name, const std::string& work_directory, std::string& file_path) const;

  // Checks for the existence of a property.
  bool checkProperty(const std::string& property_name) const;

  void treeTraversal() const;

private:

  std::unique_ptr<PropertyImpl> properties_impl_ptr_;    // boost:: properties PIMPL

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold genome auxiliary file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AuxFileProperty {

public:

  AuxFileProperty(const std::string& file_name, const std::string& aux_type) : file_name_(file_name), aux_type_(aux_type) {}
  AuxFileProperty(const AuxFileProperty&) = default;
  ~AuxFileProperty() = default;

  const std::string& fileName() const { return file_name_; }
  const std::string& auxType() const { return aux_type_; }

  bool supportedFile(const std::string& aux_type) const { return supported_types_.find(aux_type) != supported_types_.end(); }

  constexpr static const char* AUX_FILE_NAME_ = "fileName.";
  constexpr static const char* AUX_FILE_TYPE_ = "auxType.";

  // Supported auxiliary file types.
  // GFF file that contains the offsets for the Translation Start Sequences (TSS) in the Genome of Plasmodium Falciparum (3D7).
  constexpr static const char* ADJALLEY_TSS_GFF_ = "ADJALLEY_TSS_GFF";

private:

  std::string file_name_;
  std::string aux_type_;

  const std::set<std::string> supported_types_ = { ADJALLEY_TSS_GFF_ };


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// High level object extracts application specific properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeProperties {

public:

  RuntimeProperties() = default;
  ~RuntimeProperties() = default;

  bool readProperties(const std::string& properties_file);

  void setWorkDirectory(const std::string& work_directory) { work_directory_ = work_directory; }

  const std::string& workDirectory() const { return work_directory_; }

  bool getMixtureFile(std::string& mixture_file) const;

  bool getPropertiesAuxFile(std::string &aux_file) const;

  size_t getVCFPloidy() const;

  void getGenomeDBFiles(const std::string& organism,
                        std::string& fasta_file,
                        std::string& gff_file,
                        std::string& gaf_file,
                        std::string& tranlation_table) const;

  bool getGenomeAuxFiles(const std::string& organism, std::vector<AuxFileProperty>& auxfiles) const;

  bool getVCFFiles(std::vector<std::string>& vcf_files) const;

  bool getActiveGenomes(std::vector<std::string>& genome_list) const;

private:

  // Node categories.
  constexpr static const char* RUNTIME_ROOT_ = "runTime.";
  constexpr static const char* HELP_ = ".help";
  constexpr static const char* VALUE_ = ".value";
  // Runtime categories.
  constexpr static const char* VCF_LIST_ = "vcfList.";
  constexpr static const char* GENOME_LIST_ = "genomeList.";
  constexpr static const char* ACTIVE_ = "active";
  constexpr static const char* FILE_LIST_ = "fileList";
  constexpr static const char* VCF_PLOIDY_ = "vcfPloidy";
  constexpr static const char* FASTA_FILE_ = "fastaFile";
  constexpr static const char* GFF_FILE_ = "gffFile";
  constexpr static const char* GAF_FILE_ = "gafFile";
  constexpr static const char* MIXTURE_FILE_ = "mixtureFile";
  constexpr static const char* AUX_FILE_ = "auxFile";
  constexpr static const char* AUX_GENOME_FILE_ = "auxGenomeFile";
  constexpr static const char* TRANSLATION_TABLE_ = "translationTable";
  // Defaults
  constexpr static const size_t DEFAULT_PLOIDY_ = 2;

  std::string work_directory_;  // The work directory, all files are specified 'work_directory/file_name'
  PropertyTree property_tree_;   // All the option XML files.

};






}   // namespace genome
}   // namespace kellerberrin





#endif //KGL_PROPERTIES_H
