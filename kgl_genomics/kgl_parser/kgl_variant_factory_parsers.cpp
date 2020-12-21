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

kgl::DataFileParserEnum kgl::ParserSelection::getParserType(const std::string& parser_type) {

  std::string parser_upper = Utility::toupper(parser_type);

  for (auto const& [parser_type, parser_string] : implementated_parsers_) {

    if (parser_upper == Utility::toupper(parser_string)) {

      return parser_type;

    }

  }

  return DataFileParserEnum::NotImplemented;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

// Read and parse the specified ancestry file.
std::shared_ptr<kgl::DataObjectBase> kgl::ParserSelection::readPEDAncestry(std::shared_ptr<BaseFileInfo> file_info) {

  auto ped_file_info = std::dynamic_pointer_cast<PedAncestryInfo>(file_info);

  if (not ped_file_info) {

    ExecEnv::log().critical("ParserSelection::readPEDAncestry, Expected PED (ancestry .ped) file for file ident: {}", file_info->identifier());

  }

  std::shared_ptr<GenomePEDData> ped_data(std::make_shared<GenomePEDData>(ped_file_info->identifier()));

  ParsePedFile ped_parser(ped_data);

  ped_parser.readParsePEDImpl(ped_file_info->fileName());

  return ped_data;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

std::shared_ptr<kgl::DataObjectBase> kgl::ParserSelection::parseData( std::shared_ptr<const GenomeCollection> reference_genomes,
                                                                      std::shared_ptr<BaseFileInfo> file_info_ptr,
                                                                      const VariantEvidenceMap& evidence_map,
                                                                      const ContigAliasMap& contig_alias) {

  switch(ParserSelection::getParserType(file_info_ptr->parserIdent())) {

    case DataFileParserEnum::GatkMultiGenome:
      return readVCF<Pf3kVCFImpl, UnphasedPopulation>(reference_genomes, file_info_ptr, evidence_map, contig_alias);

    case DataFileParserEnum::GRChNoGenome:
      return readVCF<GrchVCFImpl, UnphasedPopulation>(reference_genomes, file_info_ptr, evidence_map, contig_alias);

    case DataFileParserEnum::MultiGenomePhased:
      return readVCF<Genome1000VCFImpl, DiploidPopulation>(reference_genomes, file_info_ptr, evidence_map, contig_alias);

    case DataFileParserEnum::MultiGenomeGnomad:
      return readVCF<GenomeGnomadVCFImpl, DiploidPopulation>(reference_genomes, file_info_ptr, evidence_map, contig_alias);

    case DataFileParserEnum::PedAncestry:
      return ParserSelection::readPEDAncestry(file_info_ptr);

    case DataFileParserEnum::NotImplemented:
      ExecEnv::log().critical("ParserSelection::parseData; Package: {}, data file ident: {}, parser: {} not defined",
                              file_info_ptr->fileName(), file_info_ptr->parserIdent());

  }

  return ParserSelection::readPEDAncestry(file_info_ptr); // Never reached.

}