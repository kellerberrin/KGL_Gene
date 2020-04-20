//
// Created by kellerberrin on 11/11/18.
//


#ifndef KGL_PROPERTIES_H
#define KGL_PROPERTIES_H

#include "kel_property_tree.h"

#include <memory>
#include <string>
#include <set>


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold genome auxiliary file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AuxFileInfo {

public:

  AuxFileInfo(const std::string& file_name, const std::string& aux_type) : file_name_(file_name), aux_type_(aux_type) {}
  AuxFileInfo(const AuxFileInfo&) = default;
  ~AuxFileInfo() = default;

  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  [[nodiscard]] const std::string& auxType() const { return aux_type_; }

  [[nodiscard]] bool supportedFile(const std::string& aux_type) const { return supported_types_.find(aux_type) != supported_types_.end(); }

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
// Object to vcf file infomation.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class VCFParserEnum{ GatkMultiGenome, GRChNoGenome, NotImplemented};
using VCFParserTypes = std::vector<std::pair<VCFParserEnum, std::string>>;

class VCFFileInfo {

public:

  VCFFileInfo(const std::string& file_name, const std::string& parser_type, const std::string& reference_genome, size_t ploidy)
      : file_name_(file_name),
        reference_genome_(reference_genome),
        ploidy_(ploidy) { parser_type_ = getParserType(parser_type); }
  VCFFileInfo(const VCFFileInfo&) = default;
  ~VCFFileInfo() = default;

  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  [[nodiscard]] const std::string& referenceGenome() const { return reference_genome_; }
  [[nodiscard]] VCFParserEnum parserType() const { return parser_type_; }
  [[nodiscard]] size_t ploidy() const { return ploidy_; }

  constexpr static const char* AUX_FILE_NAME_ = "fileName.";
  constexpr static const char* AUX_FILE_TYPE_ = "auxType.";


private:


  std::string file_name_;
  std::string reference_genome_;
  size_t ploidy_;
  VCFParserEnum parser_type_;
  const VCFParserTypes implementated_parsers_{ std::pair<VCFParserEnum, std::string>(VCFParserEnum::GatkMultiGenome, "GatkMultiGenome"),
                                               std::pair<VCFParserEnum, std::string>(VCFParserEnum::GRChNoGenome, "GRChNoGenome"),
                                               std::pair<VCFParserEnum, std::string>(VCFParserEnum::NotImplemented, "NotImplemented")};

  VCFParserEnum getParserType(const std::string& parser_type) const;

};

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

  [[nodiscard]] bool getMixtureFile(std::string& mixture_file) const;

  [[nodiscard]] bool getPropertiesAuxFile(std::string &aux_file) const;


  void getGenomeDBFiles(const std::string& organism,
                        std::string& fasta_file,
                        std::string& gff_file,
                        std::string& gaf_file,
                        std::string& tranlation_table) const;

  [[nodiscard]] bool getGenomeAuxFiles(const std::string& organism, std::vector<AuxFileInfo>& auxfiles) const;


  [[nodiscard]] bool getActiveGenomes(std::vector<std::string>& genome_list) const;

  [[nodiscard]] std::vector<VCFFileInfo> getVCFFileVector() const;

private:

  constexpr static const char* DOT_ = ".";
  // Node categories.
  constexpr static const char* RUNTIME_ROOT_ = "runTime";
  constexpr static const char* HELP_ = "help";
  constexpr static const char* ACTIVE_ = "active";
  constexpr static const char* VALUE_ = "value";
  // Runtime categories.
  constexpr static const char* VCF_LIST_ = "vcfList";
  constexpr static const char* VCF_FILE_ = "vcfFile";
  constexpr static const char* VCF_FILE_NAME_ = "vcfFileName";
  constexpr static const char* VCF_PARSER_TYPE_ = "vcfParser";
  constexpr static const char* VCF_FILE_GENOME_ =  "vcfGenome";
  constexpr static const char* VCF_FILE_PLOIDY_ =  "vcfPloidy";
  constexpr static const char* GENOME_LIST_ = "genomeList";
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


}   // end namespace



#endif //KGL_PROPERTIES_H
