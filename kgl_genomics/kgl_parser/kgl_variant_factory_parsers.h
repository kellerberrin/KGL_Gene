//
// Created by kellerberrin on 20/12/20.
//

#ifndef KGL_VARIANT_FACTORY_PARSERS_H
#define KGL_VARIANT_FACTORY_PARSERS_H

#include "kgl_runtime.h"
#include "kgl_resource_db.h"



namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Provides a common interface for data parsers.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class ParserSelection {

public:

  ParserSelection() = default;
  ~ParserSelection() = default;

  [[nodiscard]] static std::shared_ptr<DataDB> parseData(const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                                         std::shared_ptr<BaseFileInfo> file_info,
                                                         const VariantEvidenceMap& evidence_map,
                                                         const ContigAliasMap& contig_alias);

private:


  // Called when the PedAncestry parser is specified to read in an ancestry (.ped) file
  [[nodiscard]] static std::shared_ptr<DataDB> readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info,
                                                               DataSourceEnum data_source);

  // Read Pf3k Complexity of Infection information.
  [[nodiscard]] static std::shared_ptr<DataDB> readPf3kCOI( std::shared_ptr<BaseFileInfo> file_info,
                                                            DataSourceEnum data_source);

  // Read and parse a package specified VCF file.
  template<class VCFParser>
  [[nodiscard]] static std::shared_ptr<DataDB> readVCF(const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                                       std::shared_ptr<BaseFileInfo> file_info,
                                                       const VariantEvidenceMap& evidence_map,
                                                       const ContigAliasMap& contig_alias,
                                                       DataSourceEnum data_source) {

    // Get the physical file name, VCF file type etc.
    auto vcf_file_info = std::dynamic_pointer_cast<RuntimeVCFFileInfo>(file_info);

    if (not vcf_file_info) {

      ExecEnv::log().critical("ExecutePackage::readVCF; Expected VCF file for file ident: {}", file_info->identifier());

    }

    // Get the specified reference genome to validate the parsed VCF population.
    std::optional<std::shared_ptr<const GenomeReference>> ref_genome_opt = resource_ptr->getGenomes().getOptionalResource(vcf_file_info->referenceGenome());

    if (not ref_genome_opt) {

      ExecEnv::log().critical("ParserSelection::readVCF; Reference Genome {} Not Found for VCF file ident: {}",
                              vcf_file_info->referenceGenome(), vcf_file_info->identifier());

    }

    std::shared_ptr<const GenomeReference> ref_genome = ref_genome_opt.value();

    // Get the defined INFO subset defined for the VCF (can be all INFO fields).
    auto evidence_opt = evidence_map.lookupEvidence(vcf_file_info->evidenceIdent());

    if (not evidence_opt) {

      ExecEnv::log().critical("ParserSelection::readVCF; Evidence Ident {} Not Found for VCF file ident: {}",
                              vcf_file_info->evidenceIdent(), vcf_file_info->identifier());

    }

    // The variant population and VCF data source.
    std::shared_ptr<PopulationDB> vcf_population_ptr(std::make_shared<PopulationDB>(vcf_file_info->identifier(), data_source));

    // Read the VCF with the appropriate parser in a unique block so that that the parser is deleted before validation begins.
    // This prevents the parser queues stall warning from activating if the population verification is lengthy.
    {
      VCFParser reader(vcf_population_ptr, ref_genome_opt.value(), contig_alias, evidence_opt.value());
      reader.readParseVCFImpl(vcf_file_info->fileName());
    }

    // Validate the parsed VCF population against the specified reference genome.
    auto [total_variants, validated_variants] = vcf_population_ptr->validate(ref_genome);

    ExecEnv::log().info("File: {}, Total Variants: {}, Validated Variants: {} ({})",
                        vcf_population_ptr->populationId(), total_variants, validated_variants, ref_genome->genomeId());

    return vcf_population_ptr;

  }


};




} // namespace

#endif // KGL_VARIANT_FACTORY_PARSERS_H
