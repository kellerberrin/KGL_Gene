//
// Created by kellerberrin on 31/03/23.
//


#include "kgl_properties.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;


// A map of resources.
kgl::ResourceDefinitions kgl::RuntimeProperties::getRuntimeResourceDef() const {

  ResourceDefinitions resource_definitions;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + RESOURCE_LIST_;
  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().warn("RuntimeProperties::getRuntimeResources; No runtime resources specified.");
    return resource_definitions; // return empty map.

  }

  for (const auto &[tree_type, sub_tree]: property_tree_vector) {

    std::optional<ResourceParameters> resource_opt{std::nullopt};
    switch (Utility::hash(tree_type)) {

      case Utility::hash(GENOME_DATABASE_):
        resource_opt = genomeResource(sub_tree);
        break;

      case Utility::hash(ONTOLOGY_DATABASE_):
        resource_opt = ontologyDatabase(sub_tree);
        break;

      case Utility::hash(GENE_ID_DATABASE_):
        resource_opt = geneIDDatabase(sub_tree);
        break;

      case Utility::hash(GENEALOGY_ID_DATABASE_):
        resource_opt = genealogyIDDatabase(sub_tree);
        break;

      case Utility::hash(CITATION_DATABASE_):
        resource_opt = citationDatabase(sub_tree);
        break;

      case Utility::hash(ENTREZ_DATABASE_):
        resource_opt = entrezDatabase(sub_tree);
        break;

      case Utility::hash(PF7_SAMPLE_DATABASE_):
        resource_opt = Pf7SampleDatabase(sub_tree);
        break;

      case Utility::hash(PUBMED_LIT_API_):
        resource_opt = PubmedLitAPI(sub_tree);
        break;

      case Utility::hash(AUX_ID_DATABASE_):
        resource_opt = auxIDDatabase(sub_tree);
        break;

      case Utility::hash(COMMENT_):
        break;

      default:
        ExecEnv::log().warn("RuntimeProperties::getRuntimeResources, Unexpected resource: '{}'; ignored", tree_type);
        break;

    }

    if (resource_opt) {

      resource_definitions.insert({tree_type, resource_opt.value()});

    }

  }

  return resource_definitions;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::genomeResource(const PropertyTree& sub_tree) const {

  std::string genome_ident;
  if (not sub_tree.getProperty(GENOME_IDENT_, genome_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Database Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(GENOME_DATABASE_, genome_ident);

  std::string fasta_file_name;
  if (not sub_tree.getFileProperty(FASTA_FILE_, workDirectory(), fasta_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Fasta file name information.");
    return std::nullopt;

  }
  resource_parameters.setParameter(FASTA_FILE_, fasta_file_name);

  std::string gff_file_name;
  if (not sub_tree.getFileProperty(GFF_FILE_, workDirectory(), gff_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No Gff file name information.");
    return std::nullopt;

  }
  resource_parameters.setParameter(GFF_FILE_, gff_file_name);

  std::string translation_table;
  if (not sub_tree.getProperty(TRANSLATION_TABLE_, translation_table)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No DNA/Amino translation table specified.");
    return std::nullopt;

  }
  resource_parameters.setParameter(TRANSLATION_TABLE_, translation_table);

  // The Gaf file is optional.
  std::string gaf_file_name;
  if (sub_tree.getOptionalFileProperty(GAF_ANNOTATION_FILE_, workDirectory(), gaf_file_name)) {

    resource_parameters.setParameter(GAF_ANNOTATION_FILE_, gaf_file_name);

  }

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::ontologyDatabase(const PropertyTree& sub_tree) const {

  std::string ontology_ident;
  if (not sub_tree.getProperty(ONTOLOGY_IDENT_, ontology_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Ontology Database Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(ONTOLOGY_DATABASE_, ontology_ident);

  std::string annotation_file_name;
  if (not sub_tree.getFileProperty(GAF_ANNOTATION_FILE_, workDirectory(), annotation_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Ontology Annotation (GAF) file name information.");
    return std::nullopt;

  }
  resource_parameters.setParameter(GAF_ANNOTATION_FILE_, annotation_file_name);

  std::string go_graph_file_name;
  if (not sub_tree.getFileProperty(GO_ONTOLOGY_FILE_, workDirectory(), go_graph_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No Ontology GO file name information.");
    return std::nullopt;

  }
  resource_parameters.setParameter(GO_ONTOLOGY_FILE_, go_graph_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::geneIDDatabase(const PropertyTree& sub_tree) const {


  std::string nomenclature_ident;
  if (not sub_tree.getProperty(GENE_ID_IDENT_, nomenclature_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No gene Nomenclature Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(GENE_ID_DATABASE_, nomenclature_ident);

  std::string nomenclature_file_name;
  if (not sub_tree.getFileProperty(GENE_ID_FILE_, workDirectory(), nomenclature_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Gene Nomenclature file name information, ident: {}", nomenclature_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(GENE_ID_FILE_, nomenclature_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::genealogyIDDatabase(const PropertyTree& sub_tree) const {

  std::string genealogy_ident;
  if (not sub_tree.getProperty(GENEALOGY_ID_IDENT_, genealogy_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Genealogy Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(GENEALOGY_ID_DATABASE_, genealogy_ident);

  std::string genealogy_file_name;
  if (not sub_tree.getFileProperty(GENEALOGY_ID_FILE_, workDirectory(), genealogy_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Genealogy file name information, ident: {}", genealogy_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(GENEALOGY_ID_FILE_, genealogy_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::citationDatabase(const PropertyTree& sub_tree) const {

  std::string citation_ident;
  if (not sub_tree.getProperty(CITATION_IDENT_, citation_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Allele Citation Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(CITATION_DATABASE_, citation_ident);

  std::string citation_file_name;
  if (not sub_tree.getFileProperty(CITATION_FILE_, workDirectory(), citation_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Allele Citation file name information, ident: {}", citation_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(CITATION_FILE_, citation_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::entrezDatabase(const PropertyTree& sub_tree) const {

  std::string entrez_ident;
  if (not sub_tree.getProperty(ENTREZ_IDENT_, entrez_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Entrez Gene Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(ENTREZ_DATABASE_, entrez_ident);

  std::string entrez_file_name;
  if (not sub_tree.getFileProperty(ENTREZ_FILE_, workDirectory(), entrez_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Entrez Gene file name information, ident: {}", entrez_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(ENTREZ_FILE_, entrez_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::Pf7SampleDatabase(const PropertyTree& sub_tree) const {

  std::string Pf7Sample_ident;
  if (not sub_tree.getProperty(PF7_SAMPLE_IDENT_, Pf7Sample_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(PF7_SAMPLE_DATABASE_, Pf7Sample_ident);

  std::string Pf7Sample_file_name;
  if (not sub_tree.getFileProperty(PF7_SAMPLE_FILE_, workDirectory(), Pf7Sample_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample data file name information, ident: {}", Pf7Sample_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PF7_SAMPLE_FILE_, Pf7Sample_file_name);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::PubmedLitAPI(const PropertyTree& sub_tree) const {

  std::string pubmed_api_ident;
  if (not sub_tree.getProperty(PUBMED_LIT_IDENT_, pubmed_api_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pubmed API Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(PUBMED_LIT_API_, pubmed_api_ident);

  std::string publication_cache_file;  // Note that the cache file(s) may not exist.
  if (not sub_tree.getFileCreateProperty(PUBMED_PUBLICATION_CACHE_, workDirectory(), publication_cache_file)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; Cannot create Pubmed publication API Cache file, ident: {}", pubmed_api_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PUBMED_PUBLICATION_CACHE_, publication_cache_file);

  std::string citation_cache_file;  // Note that the cache file(s) may not exist.
  if (not sub_tree.getFileCreateProperty(PUBMED_CITATION_CACHE_, workDirectory(), citation_cache_file)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; Cannot create Pubmed citation API Cache file, ident: {}", pubmed_api_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(PUBMED_CITATION_CACHE_, citation_cache_file);

  return resource_parameters;

}

std::optional<kgl::ResourceParameters> kgl::RuntimeProperties::auxIDDatabase(const PropertyTree& sub_tree) const {

  std::string genome_aux_ident;
  if (not sub_tree.getProperty(AUX_ID_IDENT_, genome_aux_ident)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Aux Identifier.");
    return std::nullopt;

  }
  ResourceParameters resource_parameters(AUX_ID_IDENT_, genome_aux_ident);

  std::string genome_aux_file_name;
  if (not sub_tree.getFileProperty(AUX_ID_FILE_, workDirectory(), genome_aux_file_name)) {

    ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Aux file name information, ident: {}", genome_aux_ident);
    return std::nullopt;

  }
  resource_parameters.setParameter(AUX_ID_FILE_, genome_aux_file_name);

  return resource_parameters;

}

