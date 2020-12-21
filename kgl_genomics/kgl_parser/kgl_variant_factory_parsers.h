//
// Created by kellerberrin on 20/12/20.
//

#ifndef KGL_VARIANT_FACTORY_PARSERS_H
#define KGL_VARIANT_FACTORY_PARSERS_H

#include "kgl_runtime.h"



namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Provides a common interface for data parsers.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class DataFileParserEnum{ GatkMultiGenome, GRChNoGenome, MultiGenomePhased, MultiGenomeGnomad, PedAncestry, NotImplemented};

class ParserSelection {

public:

  ParserSelection() = default;
  ~ParserSelection() = default;


  [[nodiscard]] static std::shared_ptr<DataObjectBase> parseData( std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                  std::shared_ptr<BaseFileInfo> file_info,
                                                                  const VariantEvidenceMap& evidence_map,
                                                                  const ContigAliasMap& contig_alias);

private:


  using VCFParserTypes = std::vector<std::pair<DataFileParserEnum, std::string>>;
  inline const static VCFParserTypes implementated_parsers_{ std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::GatkMultiGenome, "GatkMultiGenome"),
                                                             std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::GRChNoGenome, "GRChNoGenome"),
                                                             std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::MultiGenomePhased, "MultiGenomePhased"),
                                                             std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::MultiGenomeGnomad, "MultiGenomeGnomad"),
                                                             std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::PedAncestry, "PedAncestry"),
                                                             std::pair<DataFileParserEnum, std::string>(DataFileParserEnum::NotImplemented, "NotImplemented")};

  [[nodiscard]] static DataFileParserEnum getParserType(const std::string& parser_type);
  // Called when the PedAncestry parser is specified to read in an ancestry (.ped) file
  [[nodiscard]] static std::shared_ptr<DataObjectBase> readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info);

  template<class VCFParser, class VCFPopulation>
  [[nodiscard]] static std::shared_ptr<DataObjectBase> readVCF( std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                std::shared_ptr<BaseFileInfo> file_info,
                                                                const VariantEvidenceMap& evidence_map,
                                                                const ContigAliasMap& contig_alias) {

    auto vcf_file_info = std::dynamic_pointer_cast<RuntimeVCFFileInfo>(file_info);

    if (not vcf_file_info) {

      ExecEnv::log().critical("ExecutePackage::readVCF; Expected VCF file for file ident: {}", file_info->identifier());

    }

    std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = reference_genomes->getOptionalGenome(vcf_file_info->referenceGenome());

    if (not ref_genome_opt) {

      ExecEnv::log().critical("ParserSelection::readVCF; Reference Genome {} Not Found for VCF file ident: {}",
                              vcf_file_info->referenceGenome(), vcf_file_info->identifier());

    }

    auto evidence_opt = evidence_map.lookupEvidence(vcf_file_info->evidenceIdent());

    if (not evidence_opt) {

      ExecEnv::log().critical("ParserSelection::readVCF; Evidence Ident {} Not Found for VCF file ident: {}",
                              vcf_file_info->evidenceIdent(), vcf_file_info->identifier());

    }

    // Read variants.
    std::shared_ptr<VCFPopulation> vcf_population_ptr(std::make_shared<VCFPopulation>(vcf_file_info->identifier()));

    VCFParser reader(vcf_population_ptr, ref_genome_opt.value(), contig_alias, evidence_opt.value());
    reader.readParseVCFImpl(vcf_file_info->fileName());

    auto [total_variants, validated_variants] = vcf_population_ptr->validate(ref_genome_opt.value());

    ExecEnv::log().info("Genome: {}, Total Variants: {}, Validated Variants: {}", vcf_population_ptr->populationId(), total_variants, validated_variants);

    return vcf_population_ptr;

  }


};




} // namespace

#endif // KGL_VARIANT_FACTORY_PARSERS_H
