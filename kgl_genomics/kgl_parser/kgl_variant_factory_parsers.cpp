//
// Created by kellerberrin on 20/12/20.
//

#include "kgl_ped_parser.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_factory_1000_impl.h"
#include "kgl_variant_factory_gnomad_impl.h"


#include "kgl_variant_factory_parsers.h"


namespace kgl = kellerberrin::genome;



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

// Read and parse the specified ancestry file.
std::shared_ptr<kgl::DataDB> kgl::ParserSelection::readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info, DataSourceEnum data_source) {

  auto ped_file_info = std::dynamic_pointer_cast<PedAncestryInfo>(file_info);

  if (not ped_file_info) {

    ExecEnv::log().critical("ParserSelection::readPEDAncestry, Expected PED (ancestry .ped) file for file ident: {}", file_info->identifier());

  }

  std::shared_ptr<GenomePEDData> ped_data(std::make_shared<GenomePEDData>(ped_file_info->identifier(), data_source));

  ParsePedFile ped_parser(ped_data);

  ped_parser.readParsePEDImpl(ped_file_info->fileName());

  return ped_data;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


std::shared_ptr<kgl::DataDB> kgl::ParserSelection::parseData(std::shared_ptr<const GenomeCollection> reference_genomes,
                                                             std::shared_ptr<BaseFileInfo> file_info_ptr,
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
      return readVCF<Pf3kVCFImpl>(reference_genomes, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::MonoGenomeUnphased:
      return readVCF<GrchVCFImpl>(reference_genomes, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidPhased:
      return readVCF<Genome1000VCFImpl>(reference_genomes, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidGnomad:
      return readVCF<GenomeGnomadVCFImpl>(reference_genomes, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::PedGenome1000:
      return ParserSelection::readPEDAncestry(file_info_ptr, data_source);

  }

  // Never reached.
  return ParserSelection::readPEDAncestry(file_info_ptr, data_source);

}