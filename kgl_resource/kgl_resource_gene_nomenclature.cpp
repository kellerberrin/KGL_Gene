//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_gene_nomenclature.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_uniprot_parser.h"
#include "kgl_ensembl_id_parser.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view NOMENCLATURE_IDENT{"nomenclatureIdent"};
  constexpr std::string_view NOMENCLATURE_FILE{"nomenclatureFile"};
}


/// Returns the gene nomenclature resource type identifier.
std::string_view kgl::GeneNomenclatureFactory::resourceType() const {

  return ResourceType::GENE_NOMENCLATURE;

}


/// Parse the gene nomenclature resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::GeneNomenclatureFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "GeneNomenclatureFactory::parseXML");
  reader.requiredString(NOMENCLATURE_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::GENE_NOMENCLATURE), ident);
  std::string file_name;
  reader.requiredFile(NOMENCLATURE_FILE, work_directory, file_name);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(NOMENCLATURE_FILE), std::move(file_name));
  return params;

}


/// Construct the gene nomenclature resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::GeneNomenclatureFactory::load(const ResourceParameters& parameters) const {

  auto file_opt = parameters.getParameter(std::string(NOMENCLATURE_FILE));
  if (not file_opt) {

    ExecEnv::log().critical("GeneNomenclatureFactory::load; nomenclature file not defined for: {}",
                            parameters.resourceIdent());

  }

  auto const& file_name = file_opt.value();
  auto const& ident = parameters.resourceIdent();

  if (ident == ResourceType::NOMENCLATURE_UNIPROT) {

    ParseUniprotId parse_uniprot;
    if (not parse_uniprot.parseUniprotFile(file_name)) {

      ExecEnv::log().critical("GeneNomenclatureFactory::load; failed to parse Uniprot file: {}", file_name);

    }
    return std::make_shared<const UniprotResource>(ident, parse_uniprot.getUniproResource());

  }

  if (ident == ResourceType::NOMENCLATURE_ENSEMBL) {

    ParseGeneIdents parse_gene_idents;
    if (not parse_gene_idents.parseIdentFile(file_name)) {

      ExecEnv::log().critical("GeneNomenclatureFactory::load; failed to parse Ensembl file: {}", file_name);

    }
    return std::make_shared<const EnsemblHGNCResource>(ident, parse_gene_idents.getSynonymVector());

  }

  ExecEnv::log().critical("GeneNomenclatureFactory::load; unknown nomenclature ident: {}", ident);

}
