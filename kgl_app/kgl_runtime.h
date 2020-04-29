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



class VCFFileInfo;
using VCFFileMap = std::map<std::string, VCFFileInfo>;
enum class VCFParserEnum{ GatkMultiGenome, GRChNoGenome, NotImplemented};

class VCFFileInfo {

public:

  VCFFileInfo(const std::string& vcf_identifier,
              const std::string& file_name,
              const std::string& parser_type,
              const std::string& reference_genome,
              size_t ploidy)
  : vcf_identifier_(vcf_identifier),
    file_name_(file_name),
    reference_genome_(reference_genome),
    ploidy_(ploidy) { parser_type_ = getParserType(parser_type); }
  VCFFileInfo(const VCFFileInfo&) = default;
  ~VCFFileInfo() = default;

  [[nodiscard]] const std::string& identifier() const { return vcf_identifier_; }
  [[nodiscard]] const std::string& fileName() const { return file_name_; }
  [[nodiscard]] const std::string& referenceGenome() const { return reference_genome_; }
  [[nodiscard]] VCFParserEnum parserType() const { return parser_type_; }
  [[nodiscard]] size_t ploidy() const { return ploidy_; }


private:

  std::string vcf_identifier_;   // A unique short string to identify this VCF file in other classes
  std::string file_name_;
  std::string reference_genome_;
  size_t ploidy_;
  VCFParserEnum parser_type_;
  using VCFParserTypes = std::vector<std::pair<VCFParserEnum, std::string>>;
  const VCFParserTypes implementated_parsers_{ std::pair<VCFParserEnum, std::string>(VCFParserEnum::GatkMultiGenome, "GatkMultiGenome"),
                                               std::pair<VCFParserEnum, std::string>(VCFParserEnum::GRChNoGenome, "GRChNoGenome"),
                                               std::pair<VCFParserEnum, std::string>(VCFParserEnum::NotImplemented, "NotImplemented")};

  VCFParserEnum getParserType(const std::string& parser_type) const;

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

  [[nodiscard]] const ContigId_t& lookupAlias(const ContigId_t& alias);
  void setAlias(const ContigId_t& alias, const ContigId_t& contig_id);

private:

  AliasMap alias_map_;

};


} // namespace

#endif //KGL_KGL_RUNTIME_H
