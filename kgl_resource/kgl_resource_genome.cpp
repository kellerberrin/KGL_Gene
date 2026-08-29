//
// Created by kellerberrin on 26/08/26.
//

#include "kgl_resource_genome.h"
#include "kgl_resource_type.h"
#include "kgl_xml_reader.h"
#include "kgl_genome_genome.h"

namespace kgl = kellerberrin::genome;

namespace {
  constexpr std::string_view GENOME_IDENT{"genomeIdent"};
  constexpr std::string_view FASTA_FILE{"fastaFile"};
  constexpr std::string_view GFF_FILE{"gffFile"};
  constexpr std::string_view TRANSLATION_TABLE{"translationTable"};
  constexpr std::string_view GAF_ANNOTATION_FILE{"gafFile"};
}


/// Returns the genome database resource type identifier.
std::string_view kgl::GenomeResourceFactory::resourceType() const {

  return ResourceType::GENOME;

}


/// Parse the genome resource XML sub-tree into ResourceParameters.
std::optional<kgl::ResourceParameters>
kgl::GenomeResourceFactory::parseXML(const PropertyTree& sub_tree, std::string_view work_directory) const {

  std::string ident;
  XmlReader reader(sub_tree, "GenomeResourceFactory::parseXML");
  reader.requiredString(GENOME_IDENT, ident);

  ResourceParameters params(std::string(ResourceType::GENOME), ident);
  std::string fasta, gff, translation, gaf;
  reader.requiredFile(FASTA_FILE, work_directory, fasta)
        .requiredFile(GFF_FILE, work_directory, gff)
        .requiredString(TRANSLATION_TABLE, translation);
  if (not reader.ok()) return std::nullopt;

  params.setParameter(std::string(FASTA_FILE), std::move(fasta));
  params.setParameter(std::string(GFF_FILE), std::move(gff));
  params.setParameter(std::string(TRANSLATION_TABLE), std::move(translation));
  if (reader.optionalFile(GAF_ANNOTATION_FILE, work_directory, gaf)) {

    params.setParameter(std::string(GAF_ANNOTATION_FILE), std::move(gaf));

  }
  return params;

}


/// Construct the genome database resource from its parameters.
std::shared_ptr<const kgl::ResourceBase>
kgl::GenomeResourceFactory::load(const ResourceParameters& parameters) const {

  auto fasta_opt = parameters.getParameter(std::string(FASTA_FILE));
  auto gff_opt = parameters.getParameter(std::string(GFF_FILE));
  auto trans_opt = parameters.getParameter(std::string(TRANSLATION_TABLE));
  if (not fasta_opt or not gff_opt or not trans_opt) {

    ExecEnv::log().critical("GenomeResourceFactory::load; missing required parameters for genome: {}",
                            parameters.resourceIdent());

  }

  std::string gaf_file;
  if (auto gaf_opt = parameters.getParameter(std::string(GAF_ANNOTATION_FILE))) {

    gaf_file = gaf_opt.value();

  }

  return GenomeReference::createGenomeDatabase(parameters.resourceIdent(),
                                                fasta_opt.value(),
                                                gff_opt.value(),
                                                gaf_file,
                                                trans_opt.value());

}