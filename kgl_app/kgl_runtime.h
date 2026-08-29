//
// Created by kellerberrin on 29/4/20.
//

#ifndef KGL_KGL_RUNTIME_H
#define KGL_KGL_RUNTIME_H

#include "kgl_genome_types.h"
#include "kgl_runtime_resource.h"

#include <limits>
#include <memory>
#include <string>
#include <vector>
#include <optional>
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

/// Holds a package identifier for an active package that will be executed at runtime.
class ActivePackage {

public:

  ActivePackage(const std::string& package_identifier) : package_identifier_(package_identifier) {}
  ~ActivePackage() = default;

  /// Returns the package identifier string.
  [[nodiscard]] const std::string& packageIdentifier() const { return package_identifier_; }

private:

  std::string package_identifier_;


};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold a Package object and associated analysis, genome database and VCF file objects (only identifiers).
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RuntimePackage;
using RuntimePackageMap = std::map<std::string, RuntimePackage>;

/// Holds a package definition with associated analysis, resource, and iterative file lists.
class RuntimePackage {

public:

  RuntimePackage( std::string package_identifier,
                  std::vector<std::string> analysis_list,
                  std::vector<std::pair<std::string, std::string>> resource_list,
                  std::vector<std::vector<std::string>> iterative_file_list)
                  : package_identifier_(std::move(package_identifier)),
                    analysis_list_(std::move(analysis_list)),
                    resource_list_(std::move(resource_list)),
                    iterative_file_list_(std::move(iterative_file_list)) {}
  RuntimePackage(const RuntimePackage&) = default;
  ~RuntimePackage() = default;

  /// Returns the package identifier string.
  [[nodiscard]] const std::string& packageIdentifier() const { return package_identifier_; }
  /// Returns the list of analysis identifiers associated with this package.
  [[nodiscard]] const std::vector<std::string>& analysisList() const { return analysis_list_; }
  /// Returns the list of resource identifier/type pairs associated with this package.
  [[nodiscard]] const std::vector<std::pair<std::string, std::string>>& resourceList() const { return resource_list_; }
  /// Returns the iterative file list (groups of file identifiers) for this package.
  [[nodiscard]] const std::vector<std::vector<std::string>>& iterativeFileList() const { return iterative_file_list_; }

private:

  std::string package_identifier_;
  std::vector<std::string> analysis_list_;
  std::vector<std::pair<std::string, std::string>> resource_list_;
  std::vector<std::vector<std::string>> iterative_file_list_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold Analysis object and the associated parameter list for each analysis.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class RuntimeAnalysis;
using RuntimeAnalysisMap = std::map<std::string, RuntimeAnalysis>;
using RuntimeParameterMap = std::vector<std::string>;

/// Holds an analysis definition and its associated parameter list.
class RuntimeAnalysis {

public:

  RuntimeAnalysis(std::string analysis_identifier, RuntimeParameterMap parameter_map)
  : analysis_identifier_(std::move(analysis_identifier)),
    parameter_map_(std::move(parameter_map)) {}
  RuntimeAnalysis(const RuntimeAnalysis&) = default;
  ~RuntimeAnalysis() = default;

  /// Returns the parameter map associated with this analysis.
  [[nodiscard]] const RuntimeParameterMap& parameterMap() const { return parameter_map_; }
  /// Returns the analysis identifier string.
  [[nodiscard]] const std::string& analysis() const { return analysis_identifier_; }

private:

  std::string analysis_identifier_;
  RuntimeParameterMap parameter_map_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Base Object to hold data file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class BaseFileInfo;
using RuntimeDataFileMap = std::map<std::string, std::shared_ptr<BaseFileInfo>>;

/// Base class holding data file information parsed from the runtime XML.
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

  /// Returns the unique file identifier string.
  [[nodiscard]] const std::string& identifier() const { return file_identifier_; }
  /// Returns the file name.
  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  /// Returns the file type string.
  [[nodiscard]] const std::string& fileType() const { return file_type_; }

private:

  const std::string file_identifier_;   // A unique short string to identify this VCF file in other classes
  const std::string file_name_;
  const std::string file_type_;

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Object to hold vcf file information.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// Holds VCF file information including the reference genome and evidence identifiers.
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

  /// Returns the reference genome identifier for this VCF file.
  [[nodiscard]] const std::string& referenceGenome() const { return reference_genome_; }
  /// Returns the evidence identifier for this VCF file.
  [[nodiscard]] const std::string& evidenceIdent() const { return evidence_ident_; }


private:

  std::string reference_genome_;
  std::string evidence_ident_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Lookup the Homosapien fasta/gff contig/chromosome  identifier using a VCF contig_ref_ptr identifier.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Enum representing the type of chromosome.
enum class ChromosomeType { AUTOSOMAL, ALLOSOME_X, ALLOSOME_Y, MITOCHONDRIA };
using AliasMap = std::unordered_map<ContigId_t, std::pair<ContigId_t, ChromosomeType>>;

/// Maps VCF contig identifiers to canonical contig identifiers and chromosome types.
class ContigAliasMap {

public:

  ContigAliasMap() = default;
  ContigAliasMap(const ContigAliasMap&) = default;
  ~ContigAliasMap() = default;

  /// Looks up the canonical contig identifier for the given alias.
  [[nodiscard]] const ContigId_t& lookupAlias(const ContigId_t& alias) const;
  /// Looks up the chromosome type for the given alias.
  [[nodiscard]] ChromosomeType lookupType(const ContigId_t& alias) const;

  /// Sets an alias mapping from a VCF contig identifier to a canonical contig identifier and chromosome type.
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

/// Indexed by the identifier <evidenceIdent>, an ordered set of INFO field IDs.
using EvidenceInfoSet = std::set<std::string>;
using EvidenceMap = std::map<std::string, EvidenceInfoSet>;

/// Holds VCF INFO evidence specifications indexed by evidence identifier.
class VariantEvidenceMap {

public:

  VariantEvidenceMap() = default;
  VariantEvidenceMap(const VariantEvidenceMap&) = default;
  ~VariantEvidenceMap() = default;

  /// Returns the underlying evidence map.
  [[nodiscard]] const EvidenceMap& getMap() const { return evidence_map_; }
  /// Looks up the set of INFO field IDs for the given evidence identifier.
  [[nodiscard]] std::optional<const EvidenceInfoSet> lookupEvidence(const std::string& evidence_ident) const;
  /// Sets the evidence specification for the given evidence identifier.
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

/// A named multimap of parameter values where parameter names need not be unique.
class ParameterMap {

public:

  ParameterMap() = default;
  ~ParameterMap() = default;

  /// Inserts a parameter with the given identifier and value.
  void insert(const std::string& ident, const std::string& value) { parameter_map_.emplace(ident, value); }
  /// Retrieves all values for the given parameter identifier.
  [[nodiscard]] std::vector<std::string> retrieve(const std::string& ident) const;

  // Returns std::nullopt if parameter is wrong size e.g. if parameter not found (0) and size specified = 1.
  // ANY_SIZE will always return a std::vector (possibly empty).
  /// Retrieves floating-point parameter values; returns std::nullopt if size mismatch.
  [[nodiscard]] std::optional<std::vector<double>> getFloat(const std::string& ident, size_t vec_size = 1) const;
  /// Retrieves string parameter values; returns std::nullopt if size mismatch.
  [[nodiscard]] std::optional<std::vector<std::string>> getString(const std::string& ident, size_t vec_size = 1) const;
  /// Retrieves integer parameter values; returns std::nullopt if size mismatch.
  [[nodiscard]] std::optional<std::vector<int64_t>> getInteger(const std::string& ident, size_t vec_size = 1) const;
  /// Retrieves size_t parameter values; returns std::nullopt if size mismatch.
  [[nodiscard]] std::optional<std::vector<size_t>> getSize(const std::string& ident, size_t vec_size = 1) const;

  /// Retrieves floating-point parameter values from a name/size pair.
  [[nodiscard]] std::optional<std::vector<double>> getFloat(const std::pair<std::string, size_t> &field) const { return getFloat(field.first, field.second); }
  /// Retrieves string parameter values from a name/size pair.
  [[nodiscard]] std::optional<std::vector<std::string>> getString(const std::pair<std::string, size_t> &field) const { return getString(field.first, field.second); }
  /// Retrieves integer parameter values from a name/size pair.
  [[nodiscard]] std::optional<std::vector<int64_t>> getInteger(const std::pair<std::string, size_t> &field) const { return getInteger(field.first, field.second); }
  /// Retrieves size_t parameter values from a name/size pair.
  [[nodiscard]] std::optional<std::vector<size_t>> getSize(const std::pair<std::string, size_t> &field) const { return getSize(field.first, field.second); }
  /// Retrieves a boolean parameter value.
  [[nodiscard]] std::optional<bool> getBool(const std::string& ident) const;

  /// Constant representing any vector size; bypasses size validation.
  constexpr static const size_t ANY_SIZE{std::numeric_limits<size_t>::max()};

private:

  // Shared implementation for typed accessors (getFloat, getInteger, getSize).
  template <typename T, typename ConvertFn>
  [[nodiscard]] std::optional<std::vector<T>> convertParameter(const std::string& ident, size_t vec_size, ConvertFn convert) const {

    std::vector<T> value_vector;
    std::vector<std::string> string_vector = retrieve(ident);

    if (string_vector.size() != vec_size and vec_size != ANY_SIZE) {
      return std::nullopt;
    }

    for (auto const& str_value : string_vector) {
      try {
        value_vector.push_back(convert(str_value));
      } catch(...) {
        ExecEnv::log().error("ParameterMap; parameter ident: {}, has invalid value: {}", ident, str_value);
        return std::nullopt;
      }
    }

    return value_vector;

  }

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

/// A unique named map of parameter vectors supplied to analysis objects.
class ActiveParameterList {

public:

  ActiveParameterList() = default;
  ~ActiveParameterList() = default;

  /// Returns the underlying parameter list map.
  [[nodiscard]] const ParameterListMap& getMap() const { return active_parameter_vectors_; }

  /// Adds a named parameter vector to the active parameter list.
  bool addNamedParameterVector(const NamedParameterVector& named_vector);
  /// Creates a sub list of named parameters based on the active parameter identifiers specified by an analysis.
  [[nodiscard]] ActiveParameterList createParameterList(const std::vector<std::string>& active_idents) const;

private:

  ParameterListMap active_parameter_vectors_;

};



} // namespace

#endif //KGL_KGL_RUNTIME_H