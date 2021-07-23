//
// Created by kellerberrin on 20/12/20.
//

#include "kgl_Hsgenealogy_parser.h"
#include "kgl_Pf3k_COI.h"
#include "kgl_variant_factory_pf_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_factory_1000_impl.h"
#include "kgl_variant_factory_gnomad_impl.h"
#include "kgl_variant_parse_json.h"


#include "kgl_variant_factory_parsers.h"


namespace kgl = kellerberrin::genome;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


std::shared_ptr<kgl::DataDB> kgl::ParserSelection::parseData(const std::shared_ptr<const AnalysisResources>& resource_ptr,
                                                             const std::shared_ptr<const BaseFileInfo>& file_info_ptr,
                                                             const VariantEvidenceMap& evidence_map,
                                                             const ContigAliasMap& contig_alias) {

  auto file_characteristic = DataDB::findCharacteristic(file_info_ptr->fileType());

  if (not file_characteristic) {

    // The file type is not defined in code.
    ExecEnv::log().critical("ParserSelection::parseData; Data file ident: {}, file type: {} not defined, cannot proceed.",
                            file_info_ptr->fileName(), file_info_ptr->fileType());

  }

  auto parser_type = file_characteristic.value().parser_type;
  auto data_source = file_characteristic.value().data_source;

  switch(parser_type) {

    case ParserTypeEnum::DiploidFalciparum:
      return readVCF<PfVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::MonoGenomeUnphased:
      return readVCF<GrchVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::MonoDBSNPUnphased:
      return readVCF<SNPdbVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidPhased:
      return readVCF<Genome1000VCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidGnomad:
      return readVCF<GenomeGnomadVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::MonoJSONdbSNPUnphased:
      return readJSONdbSNP(file_info_ptr, data_source);

    default:
      ExecEnv::log().critical("ParserSelection::parseData; Unknown data file: {} specified - unrecoverable", file_info_ptr->fileName());
      return readVCF<GenomeGnomadVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source); // never reached.

  }

}


[[nodiscard]] std::shared_ptr<kgl::DataDB> kgl::ParserSelection::readJSONdbSNP( const std::shared_ptr<const BaseFileInfo>& file_info,
                                                                                DataSourceEnum data_source) {

  // The variant population and VCF data source.
  std::shared_ptr<PopulationDB> vcf_population_ptr(std::make_shared<PopulationDB>(file_info->identifier(), data_source));

  // Read the VCF with the appropriate parser in a unique block so that that the parser is deleted before validation begins.
  // This prevents the parser queues stall warning from activating if the population verification is lengthy.
  {
    JSONInfoParser reader;
    if (reader.commenceJSONIO(file_info->fileName())) {

      size_t json_size = reader.dequeueJSONline();
      ExecEnv::log().info("JSON File: {} has size: {}", file_info->fileName(), json_size);

    }

  }

  return vcf_population_ptr;   // return the empty population for now.

}

