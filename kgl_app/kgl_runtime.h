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
#include <unordered_map>
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

  ActivePackage(const std::string& package_identifier) : package_identifier_(package_identifier) {}
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
enum class RuntimeResourceType { GENOME_DATABASE, ONTOLOGY_DATABASE, GENE_NOMENCLATURE };


class RuntimePackage {

public:

  RuntimePackage( std::string package_identifier,
                  std::vector<std::string> analysis_list,
                  std::vector<std::pair<std::string, RuntimeResourceType>> resource_database_list,
                  std::vector<std::vector<std::string>> iterative_file_list)
                  : package_identifier_(std::move(package_identifier)),
                    analysis_list_(std::move(analysis_list)),
                    resource_database_list_(std::move(resource_database_list)),
                    iterative_file_list_(std::move(iterative_file_list)) {}
  RuntimePackage(const RuntimePackage&) = default;
  ~RuntimePackage() = default;


  [[nodiscard]] const std::string& packageIdentifier() const { return package_identifier_; }
  [[nodiscard]] const std::vector<std::string>& analysisList() const { return analysis_list_; }
  [[nodiscard]] const std::vector<std::pair<std::string, RuntimeResourceType>>& resourceDatabaseList() const { return resource_database_list_; }
  [[nodiscard]] const std::vector<std::vector<std::string>>& iterativeFileList() const { return iterative_file_list_; }

private:

  std::string package_identifier_;
  std::vector<std::string> analysis_list_;
  std::vector<std::pair<std::string, RuntimeResourceType>> resource_database_list_;
  std::vector<std::vector<std::string>> iterative_file_list_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Analysis object and the associated parameter list for each analysis.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeAnalysis;
using RuntimeAnalysisMap = std::map<std::string, RuntimeAnalysis>;
using RuntimeParameterMap = std::vector<std::string>;

class RuntimeAnalysis {

public:

  RuntimeAnalysis(std::string analysis_identifier, RuntimeParameterMap parameter_map)
  : analysis_identifier_(std::move(analysis_identifier)),
    parameter_map_(std::move(parameter_map)) {}
  RuntimeAnalysis(const RuntimeAnalysis&) = default;
  ~RuntimeAnalysis() = default;

  [[nodiscard]] const RuntimeParameterMap& parameterMap() const { return parameter_map_; }
  [[nodiscard]] const std::string& analysis() const { return analysis_identifier_; }

private:

  std::string analysis_identifier_;
  RuntimeParameterMap parameter_map_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Object for resources.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RuntimeResource;
using RuntimeResourceMap = std::map<std::string, std::shared_ptr<const RuntimeResource>>;

class RuntimeResource {

public:

  RuntimeResource() = default;
  virtual ~RuntimeResource() = default;

  [[nodiscard]] virtual RuntimeResourceType resourceType() const = 0;

private:

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Genome Database file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class RuntimeGenomeResource : public RuntimeResource {

public:

  RuntimeGenomeResource(std::string genome_identifier,
                        std::string fasta_file_name,
                        std::string gff_file_name,
                        std::string translation_table)
  : genome_identifier_(std::move(genome_identifier)),
    fasta_file_name_(std::move(fasta_file_name)),
    gff_file_name_(std::move(gff_file_name)),
    translation_table_(std::move(translation_table)) {}
  RuntimeGenomeResource() = delete;
  RuntimeGenomeResource(const RuntimeGenomeResource&) = default;
  ~RuntimeGenomeResource() override = default;

  [[nodiscard]] RuntimeResourceType resourceType() const override { return RuntimeResourceType::GENOME_DATABASE; }

  [[nodiscard]] const std::string& genomeIdentifier() const { return genome_identifier_; }
  [[nodiscard]] const std::string& fastaFileName() const { return fasta_file_name_; }
  [[nodiscard]] const std::string& gffFileName() const { return gff_file_name_; }
  [[nodiscard]] const std::string& translationTable() const { return translation_table_; }
  // Gaf file and ID file are optional.
  [[nodiscard]] const std::string& gafFileName () const { return gaf_file_name_; }
  [[nodiscard]] const std::string& idFileName () const { return id_file_name_; }

  void setGafFileName(const std::string& gaf_file_name) { gaf_file_name_ = gaf_file_name; }
  void setIdFileName(const std::string& id_file_name) { id_file_name_ = id_file_name; }

private:

  std::string genome_identifier_;   // A unique short string to identify this genome database
  std::string fasta_file_name_;
  std::string gff_file_name_;
  std::string translation_table_;
  std::string gaf_file_name_;
  std::string id_file_name_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Genome Database file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeOntologyResource : public RuntimeResource {

public:

  RuntimeOntologyResource(std::string ontology_identifier,
                          std::string annotation_file_name,
                          std::string go_graph_file_name)
      : ontology_identifier_(std::move(ontology_identifier)),
        annotation_file_name_(std::move(annotation_file_name)),
        go_graph_file_name_(std::move(go_graph_file_name)) {}
  RuntimeOntologyResource() = delete;
  RuntimeOntologyResource(const RuntimeOntologyResource&) = default;
  ~RuntimeOntologyResource() override = default;

  [[nodiscard]] RuntimeResourceType resourceType() const override { return RuntimeResourceType::ONTOLOGY_DATABASE; }

  [[nodiscard]] const std::string& ontologyIdentifier() const { return ontology_identifier_; }
  [[nodiscard]] const std::string& annotationFileName() const { return annotation_file_name_; }
  [[nodiscard]] const std::string& goGraphFileName() const { return go_graph_file_name_; }

private:

  std::string ontology_identifier_;   // A unique short string to identify this annotation database
  std::string annotation_file_name_;
  std::string go_graph_file_name_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Object to hold data file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BaseFileInfo;
using RuntimeDataFileMap = std::map<std::string, std::shared_ptr<BaseFileInfo>>;


class BaseFileInfo {

public:

  BaseFileInfo(std::string identifier,
               std::string file_name,
               std::string file_type)
      : file_identifier_(std::move(identifier)),
        file_name_(std::move(file_name)),
        file_type_(std::move(file_type)) {}
  BaseFileInfo(const BaseFileInfo&) = default;
  virtual ~BaseFileInfo() = default;

  [[nodiscard]] const std::string& identifier() const { return file_identifier_; }
  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  [[nodiscard]] const std::string& fileType() const { return file_type_; }

private:

  const std::string file_identifier_;   // A unique short string to identify this VCF file in other classes
  const std::string file_name_;
  const std::string file_type_;

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold vcf file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




class RuntimeVCFFileInfo : public BaseFileInfo {

public:

  RuntimeVCFFileInfo(const std::string& identifier,
                     const std::string& file_name,
                     const std::string& file_type,
                     const std::string& reference_genome,
                     const std::string& evidence_ident)
  : BaseFileInfo(identifier, file_name, file_type),
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
// Lookup the Homosapien fasta/gff contig/chromosome  identifier using a VCF contig identifier.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class ChromosomeType { AUTOSOMAL, ALLOSOME_X, ALLOSOME_Y, MITOCHONDRIA };
using AliasMap = std::unordered_map<ContigId_t, std::pair<ContigId_t, ChromosomeType>>;

class ContigAliasMap {

public:

  ContigAliasMap() = default;
  ContigAliasMap(const ContigAliasMap&) = default;
  ~ContigAliasMap() = default;

  [[nodiscard]] const ContigId_t& lookupAlias(const ContigId_t& alias) const;
  [[nodiscard]] ChromosomeType lookupType(const ContigId_t& alias) const;

  void setAlias(const ContigId_t& alias, const ContigId_t& contig_id, const std::string& chromosome_type);

private:

  // See the alias config files for these text constants.
  constexpr static const char* AUTOSOME_ = "autosome";
  constexpr static const char* ALLOSOME_X_ = "allosomeX";
  constexpr static const char* ALLOSOME_Y_ = "allosomeY";
  constexpr static const char* MITOCHONDRIA_ = "mitochondria";

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
//
// The parameter objects, supplies named vectors of vectors of named parameters to analysis objects.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A ParameterMap is named multimap of parameter values.
// Parameter names need not be unique.

class ParameterMap {

public:

  ParameterMap() = default;
  ~ParameterMap() = default;

  void insert(const std::string& ident, const std::string& value) { parameter_map_.emplace(ident, value); }
  [[nodiscard]] std::vector<std::string> retrieve(const std::string& ident) const;

  // Returns std::nullopt if parameter is wrong size e.g. if parameter not found (0) and size specified = 1.
  // ANY_SIZE will always return a std::vector (possibly empty).
  [[nodiscard]] std::optional<std::vector<double>> getFloat(const std::string& ident, size_t vec_size = 1) const;
  [[nodiscard]] std::optional<std::vector<std::string>> getString(const std::string& ident, size_t vec_size = 1) const;
  [[nodiscard]] std::optional<std::vector<int64_t>> getInteger(const std::string& ident, size_t vec_size = 1) const;
  [[nodiscard]] std::optional<std::vector<size_t>> getSize(const std::string& ident, size_t vec_size = 1) const;

  [[nodiscard]] std::optional<std::vector<double>> getFloat(const std::pair<std::string, size_t> &field) const { return getFloat(field.first, field.second); }
  [[nodiscard]] std::optional<std::vector<std::string>> getString(const std::pair<std::string, size_t> &field) const { return getString(field.first, field.second); }
  [[nodiscard]] std::optional<std::vector<int64_t>> getInteger(const std::pair<std::string, size_t> &field) const { return getInteger(field.first, field.second); }
  [[nodiscard]] std::optional<std::vector<size_t>> getSize(const std::pair<std::string, size_t> &field) const { return getSize(field.first, field.second); }
  [[nodiscard]] std::optional<bool> getBool(const std::string& ident) const;

  constexpr static const size_t ANY_SIZE{99999999999};

private:

  std::multimap<std::string, std::string> parameter_map_;

};

// A ParameterVector is a vector of parameter maps.
// Typically a vector size of more than one will execute the underlying analysis
// multiple times with different paramter combinations.
using ParameterVector = std::vector<ParameterMap>;

// A NamedParameterVector is augmented with a unique name.
using NamedParameterVector = std::pair<std::string, ParameterVector>;

// An ActiveParameterList is unique named map of parameter vectors.
// Typically an analysis specifies one or more active NamedParameterVectors.
// Which are supplied to the analysis code as a ActiveParameterList.
// If an analysis specifies an active NamedParameterVector that does not exist
// the application terminates with an error message.
// The analysis may specify zero or more active NamedParameterVectors by name.
// The names of inactive NamedParameterVectors are conveniently parked next to the active
// NamedParameterVectors for ease of editing when deciding which parameter combination
// to run with the analysis. Note that a NamedParameterVector is usually used to
// execute the analysis multiple times.
// The ActiveParameter object is typically interpreted by the package implementation code.

using ParameterListMap = std::map<const std::string, const NamedParameterVector>;

class ActiveParameterList {

public:

  ActiveParameterList() = default;
  ~ActiveParameterList() = default;

  [[nodiscard]] const ParameterListMap& getMap() const { return active_parameter_vectors_; }

  bool addNamedParameterVector(const NamedParameterVector& named_vector);
  // Create a sub list of named parameters based on the active parameter identifiers specified by an analysis.
  [[nodiscard]] ActiveParameterList createParameterList(const std::vector<std::string>& active_idents) const;

private:

  ParameterListMap active_parameter_vectors_;

};




} // namespace

#endif //KGL_KGL_RUNTIME_H
