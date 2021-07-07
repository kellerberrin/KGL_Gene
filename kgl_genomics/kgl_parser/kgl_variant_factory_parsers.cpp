//
// Created by kellerberrin on 20/12/20.
//

#include "kgl_Hsgenealogy_parser.h"
#include "kgl_Pf3k_COI.h"
#include "kgl_variant_factory_pf_impl.h"
#include "kgl_variant_factory_grch_impl.h"
#include "kgl_variant_factory_1000_impl.h"
#include "kgl_variant_factory_gnomad_impl.h"


#include "kgl_variant_factory_parsers.h"


namespace kgl = kellerberrin::genome;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//


std::shared_ptr<kgl::DataDB> kgl::ParserSelection::parseData(const std::shared_ptr<const AnalysisResources>& resource_ptr,
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
      return readVCF<PfVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::MonoGenomeUnphased:
      return readVCF<GrchVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidPhased:
      return readVCF<Genome1000VCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    case ParserTypeEnum::DiploidGnomad:
      return readVCF<GenomeGnomadVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source);

    default:
      ExecEnv::log().critical("ParserSelection::parseData; Unknown data file: {} specified - unrecoverable", file_info_ptr->fileName());
      return readVCF<GenomeGnomadVCFImpl>(resource_ptr, file_info_ptr, evidence_map, contig_alias, data_source); // never reached.

  }

}