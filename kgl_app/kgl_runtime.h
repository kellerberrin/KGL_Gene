//
// Created by kellerberrin on 29/4/20.
//

#ifndef KGL_KGL_RUNTIME_H
#define KGL_KGL_RUNTIME_H

#include "kgl_genome_types.h"

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>

namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// These are simple intermediate objects to hold data parsed directly from the "runtime.xml" configuration file.
// These objects are passed onto the ExecPackageList object for final parsing and semantic assembly.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A simple class to hold a list of packages that will be executed at runtime.
// Note, all defined packages are loaded and parsed for correctness, but only active packages are executed.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ActivePackage;
using ActivePackageVector = std::vector<ActivePackage>;

class ActivePackage {

public:

  ActivePackage(std::string package_identifier) : package_identifier_(std::move(package_identifier)) {}
  ~ActivePackage() = default;

  [[nodiscard]] const std::string& packageIdentifier() const { return package_identifier_; }

private:

  std::string package_identifier_;


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold a Package object and associated analysis, genome database and VCF file objects (only identifiers).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RuntimePackage;
using RuntimePackageMap = std::map<std::string, RuntimePackage>;

class RuntimePackage {

public:

  RuntimePackage( const std::string& package_identifier,
                  const std::vector<std::string>& analysis_list,
                  const std::vector<std::string>& genome_database_list,
                  const std::vector<std::vector<std::string>>& iterative_file_list)
                  : package_identifier_(package_identifier),
                  analysis_list_(analysis_list),
                  genome_database_list_(genome_database_list),
                  iterative_file_list_(iterative_file_list) {}
  RuntimePackage(const RuntimePackage&) = default;
  ~RuntimePackage() = default;


  [[nodiscard]] const std::string& packageIdentifier() const { return package_identifier_; }
  [[nodiscard]] const std::vector<std::string>& analysisList() const { return analysis_list_; }
  [[nodiscard]] const std::vector<std::string>& genomeDatabaseList() const { return genome_database_list_; }
  [[nodiscard]] const std::vector<std::vector<std::string>>& iterativeFileList() const { return iterative_file_list_; }

private:

  std::string package_identifier_;
  std::vector<std::string> analysis_list_;
  std::vector<std::string> genome_database_list_;
  std::vector<std::vector<std::string>> iterative_file_list_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Analysis object and the associated parameter list for each analysis.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeAnalysis;
using RuntimeAnalysisMap = std::map<std::string, RuntimeAnalysis>;
using RuntimeParameterMap = std::map<std::string, std::string>;

class RuntimeAnalysis {

public:

  RuntimeAnalysis(const std::string& analysis_identifier, const RuntimeParameterMap& parameter_map)
  : analysis_identifier_(analysis_identifier),
    parameter_map_(parameter_map) {}
  RuntimeAnalysis(const RuntimeAnalysis&) = default;
  ~RuntimeAnalysis() = default;

  [[nodiscard]] const RuntimeParameterMap& parameterMap() const { return parameter_map_; }


private:

  std::string analysis_identifier_;
  RuntimeParameterMap parameter_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Genome Database file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeGenomeProperty;
using RuntimeGenomeDatabaseMap = std::map<std::string, RuntimeGenomeProperty>;

class RuntimeGenomeProperty {

public:

  RuntimeGenomeProperty(const std::string& genome_identifier,
                        const std::string& fasta_file_name,
                        const std::string& gff_file_name,
                        const std::string& translation_table)
  : genome_identifier_(genome_identifier),
    fasta_file_name_(fasta_file_name),
    gff_file_name_(gff_file_name),
    translation_table_(translation_table) { gaf_file_name_.clear(); }
  RuntimeGenomeProperty(const RuntimeGenomeProperty&) = default;
  ~RuntimeGenomeProperty() = default;

  [[nodiscard]] const std::string& genomeIdentifier() const { return genome_identifier_; }
  [[nodiscard]] const std::string& fastaFileName() const { return fasta_file_name_; }
  [[nodiscard]] const std::string& gffFileName() const { return gff_file_name_; }
  [[nodiscard]] const std::string& translationTable() const { return translation_table_; }
  [[nodiscard]] const std::string& gafFileName () const { return gaf_file_name_; }

  void setGafFileName(const std::string& gaf_file_name) { gaf_file_name_ = gaf_file_name; }

private:

  std::string genome_identifier_;   // A unique short string to identify this genome database
  std::string fasta_file_name_;
  std::string gff_file_name_;
  std::string translation_table_;
  std::string gaf_file_name_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Object to hold data file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class DataFileParserEnum{ GatkMultiGenome, GRChNoGenome, MultiGenomePhased, PedAncestry, NotImplemented};

class BaseFileInfo;
using RuntimeDataFileMap = std::map<std::string, std::shared_ptr<BaseFileInfo>>;


class BaseFileInfo {

public:

  BaseFileInfo(const std::string& identifier,
               const std::string& file_name,
               const std::string& parser_type)
      : file_identifier_(identifier),
        file_name_(file_name) { parser_type_ = getParserType(parser_type); }
  BaseFileInfo(const BaseFileInfo&) = default;
  virtual ~BaseFileInfo() = default;

  [[nodiscard]] const std::string& identifier() const { return file_identifier_; }
  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  [[nodiscard]] DataFileParserEnum parserType() const { return parser_type_; }


private:

  std::string file_identifier_;   // A unique short string to identify this VCF file in other classes
  std::string file_name_;
  DataFileParserEnum parser_type_;
  using VCFParserTypes = std::vector<std::pair<DataFileParserEnum, std::string>>;
  const VCFParserTypes implementated_parsers_{ std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::GatkMultiGenome, "GatkMultiGenome"),
                                               std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::GRChNoGenome, "GRChNoGenome"),
                                               std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::MultiGenomePhased, "MultiGenomePhased"),
                                               std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::PedAncestry, "PedAncestry"),
                                               std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::NotImplemented, "NotImplemented")};

  DataFileParserEnum getParserType(const std::string& parser_type) const;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold .ped ancestry file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class PedAncestryInfo : public BaseFileInfo {

public:

  PedAncestryInfo(const std::string& identifier,
                  const std::string& file_name,
                  const std::string& parser_type) : BaseFileInfo(identifier, file_name, parser_type) {}
  PedAncestryInfo(const PedAncestryInfo&) = default;
  ~PedAncestryInfo() override = default;


private:


};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold vcf file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




class RuntimeVCFFileInfo : public BaseFileInfo {

public:

  RuntimeVCFFileInfo(const std::string& identifier,
                     const std::string& file_name,
                     const std::string& parser_type,
                     const std::string& reference_genome,
                     const std::string& evidence_ident)
  : BaseFileInfo(identifier, file_name, parser_type),
    reference_genome_(reference_genome),
    evidence_ident_(evidence_ident) {}
  RuntimeVCFFileInfo(const RuntimeVCFFileInfo&) = default;
  ~RuntimeVCFFileInfo() override = default;

  [[nodiscard]] const std::string& referenceGenome() const { return reference_genome_; }
  [[nodiscard]] const std::string& evidenceIdent() const { return evidence_ident_; }


private:

  std::string reference_genome_;
  std::string evidence_ident_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lookup the fasta/gff contig/chromosome  identifier using a VCF contig identifier.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using AliasMap = std::map<ContigId_t, ContigId_t>;

class ContigAliasMap {

public:

  ContigAliasMap() = default;
  ContigAliasMap(const ContigAliasMap&) = default;
  ~ContigAliasMap() = default;

  [[nodiscard]] const AliasMap& getMap() const { return alias_map_; }
  [[nodiscard]] const ContigId_t& lookupAlias(const ContigId_t& alias) const;
  void setAlias(const ContigId_t& alias, const ContigId_t& contig_id);

private:

  AliasMap alias_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object the hold the VCF INFO evidence specification.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Indexed by the identifier <evidenceIdent>, an ordered set of INFO field IDs.
using EvidenceInfoSet = std::set<std::string>;
using EvidenceMap = std::map<std::string, EvidenceInfoSet>;

class VariantEvidenceMap {

public:

  VariantEvidenceMap() = default;
  VariantEvidenceMap(const VariantEvidenceMap&) = default;
  ~VariantEvidenceMap() = default;

  [[nodiscard]] const EvidenceMap& getMap() const { return evidence_map_; }
  [[nodiscard]] std::optional<const EvidenceInfoSet> lookupEvidence(const std::string& evidence_ident) const;
  void setEvidence(const std::string& evidence_ident, const std::set<std::string>& info_ids);

private:

  EvidenceMap evidence_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold genome auxiliary file information.
// Legacy code.
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





} // namespace

#endif //KGL_KGL_RUNTIME_H
