//
// Created by kellerberrin on 20/12/20.
//

#include "kgl_hsgenealogy_parser.h"
#include "kgl_pf3k_coi.h"
#include "kgl_variant_factory_pf_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_factory_1000_impl.h"
#include "kgl_variant_factory_gnomad_impl.h"
#include "kgl_json_parser.h"
#include "kgl_bio_pmid_parser.h"


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

    case ParserTypeEnum::FilenameOnly:
      return std::make_shared<FilenameDataDB>(data_source, file_info_ptr->fileName()); // Just return file type and file name.

    case ParserTypeEnum::ParseBioPMID:
      return readBioPMID(file_info_ptr, data_source);

    default:
      ExecEnv::log().critical("ParserSelection::parseData; Unknown data file: {} specified - unrecoverable", file_info_ptr->fileName());
      return readVCF<GenomeGnomadVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source); // never reached.

  }

}


[[nodiscard]] std::shared_ptr<kgl::DataDB> kgl::ParserSelection::readJSONdbSNP( const std::shared_ptr<const BaseFileInfo>& file_info,
                                                                                DataSourceEnum data_source) {

  // An rsid indexed map of PMID citation identifiers.
  std::shared_ptr<DBCitation> db_citation_ptr(std::make_shared<DBCitation>(data_source, file_info->fileName()));

  if (not JSONInfoParser().parseFile(file_info->fileName(), db_citation_ptr)) {

    ExecEnv::log().critical("ParserSelection::readJSONdbSNP; Error parsing Json file: {} - unrecoverable", file_info->fileName());

  }

  return db_citation_ptr;   // return the citation DB object.

}


[[nodiscard]] std::shared_ptr<kgl::DataDB> kgl::ParserSelection::readBioPMID( const std::shared_ptr<const BaseFileInfo>& file_info,
                                                                              DataSourceEnum /* data_source */) {

  // Parses a file containing Pubmed PMIDs indexed by Disease (all) and Entrez Gene ids (all species)
  ParseBioPMID bio_pmid_parser;
  if (not bio_pmid_parser.parseBioPMIDRecords(file_info->fileName())) {

    ExecEnv::log().critical("ParserSelection::readBioPMID; Error parsing Bio concepts PMID file: {} - unrecoverable", file_info->fileName());

  }
  std::shared_ptr<BioPMIDFileData> bio_pmid_ptr(std::make_shared<BioPMIDFileData>(file_info->fileName(), bio_pmid_parser.moveDiseasePMIDMap(), bio_pmid_parser.moveEntrezPMIDMap()));

  return bio_pmid_ptr;   // return the Bio PMID object.

}

