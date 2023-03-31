//
// Created by kellerberrin on 31/03/23.
//

#include "kgl_properties.h"
#include "kel_utility.h"


namespace kgl = kellerberrin::genome;



// A map of resources.
kgl::RuntimeResourceMap kgl::RuntimeProperties::getRuntimeResources() const {

  RuntimeResourceMap resource_map;

  std::string key = std::string(RUNTIME_ROOT_) + std::string(DOT_) + RESOURCE_LIST_;

  std::vector<SubPropertyTree> property_tree_vector;
  if (not property_tree_.getPropertyTreeVector(key, property_tree_vector)) {

    ExecEnv::log().warn("RuntimeProperties::getRuntimeResources; No runtime resources specified.");
    return resource_map; // return empty map.

  }

  for (const auto& [tree_type, sub_tree] : property_tree_vector) {

    // Only process Genome records.
    if (tree_type == GENOME_DATABASE_) {

      std::string genome_ident;
      if (not sub_tree.getProperty( GENOME_IDENT_, genome_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Database Identifier.");
        continue;

      }

      std::string fasta_file_name;
      if (not sub_tree.getFileProperty( FASTA_FILE_, workDirectory(), fasta_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Fasta file name information.");
        continue;

      }

      std::string gff_file_name;
      if (not sub_tree.getFileProperty( GFF_FILE_, workDirectory(), gff_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No Gff file name information.");
        continue;

      }

      std::string translation_table;
      if (not sub_tree.getProperty( TRANSLATION_TABLE_, translation_table)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No DNA/Amino translation table specified.");
        continue;

      }

      std::string gaf_file_name;
      if (not sub_tree.getOptionalFileProperty(GAF_ANNOTATION_FILE_, workDirectory(), gaf_file_name)) {

        gaf_file_name.clear();

      }


      RuntimeGenomeResource new_genome_db(genome_ident, fasta_file_name, gff_file_name, translation_table);

      new_genome_db.setGafFileName(gaf_file_name);

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeGenomeResource>(new_genome_db);

      auto const [iter, result] = resource_map.try_emplace(genome_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Genome Database ident: {} to map (duplicate)", genome_ident);

      }

    } else if (tree_type == ONTOLOGY_DATABASE_) {

      std::string ontology_ident;
      if (not sub_tree.getProperty(ONTOLOGY_IDENT_, ontology_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Ontology Database Identifier.");
        continue;

      }

      std::string annotation_file_name;
      if (not sub_tree.getFileProperty(GAF_ANNOTATION_FILE_, workDirectory(), annotation_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Ontology Annotation (GAF) file name information.");
        continue;

      }

      std::string go_graph_file_name;
      if (not sub_tree.getFileProperty(GO_ONTOLOGY_FILE_, workDirectory(), go_graph_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, No Ontology GO file name information.");
        continue;

      }

      RuntimeOntologyResource new_ontology_db(ontology_ident, annotation_file_name, go_graph_file_name);

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeOntologyResource>(new_ontology_db);

      auto const [iter, result] = resource_map.try_emplace(ontology_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Genome Database ident: {} to map (duplicate)", ontology_ident);

      }

    } else if (tree_type == GENE_ID_DATABASE_) {

      std::string nomenclature_ident;
      if (not sub_tree.getProperty(GENE_ID_IDENT_, nomenclature_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No gene Nomenclature Identifier.");
        continue;

      }

      std::string nomenclature_file_name;
      if (not sub_tree.getFileProperty(GENE_ID_FILE_, workDirectory(), nomenclature_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Gene Nomenclature file name information, ident: {}", nomenclature_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeNomenclatureResource>(nomenclature_ident,
                                                                                                                nomenclature_file_name);

      auto const [iter, result] = resource_map.try_emplace(nomenclature_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Gene Nomenclature ident: {} to map (duplicate)", nomenclature_ident);

      }

    } else if (tree_type == GENEALOGY_ID_DATABASE_) {

      std::string genealogy_ident;
      if (not sub_tree.getProperty(GENEALOGY_ID_IDENT_, genealogy_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Genealogy Identifier.");
        continue;

      }

      std::string genealogy_file_name;
      if (not sub_tree.getFileProperty(GENEALOGY_ID_FILE_, workDirectory(), genealogy_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Genealogy file name information, ident: {}", genealogy_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeGenealogyResource>(genealogy_ident, genealogy_file_name);

      auto const [iter, result] = resource_map.try_emplace(genealogy_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Genome Genealogy ident: {} to map (duplicate)", genealogy_ident);

      }

    } else if (tree_type == CITATION_DATABASE_) {

      std::string citation_ident;
      if (not sub_tree.getProperty(CITATION_IDENT_, citation_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Allele Citation Identifier.");
        continue;

      }

      std::string citation_file_name;
      if (not sub_tree.getFileProperty(CITATION_FILE_, workDirectory(), citation_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Allele Citation file name information, ident: {}", citation_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeCitationResource>(citation_ident, citation_file_name);

      auto const [iter, result] = resource_map.try_emplace(citation_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Allele Citation ident: {} to map (duplicate)", citation_ident);

      }

    } else if (tree_type == ENTREZ_DATABASE_) {

      std::string entrez_ident;
      if (not sub_tree.getProperty(ENTREZ_IDENT_, entrez_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Entrez Gene Identifier.");
        continue;

      }

      std::string entrez_file_name;
      if (not sub_tree.getFileProperty(ENTREZ_FILE_, workDirectory(), entrez_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Entrez Gene file name information, ident: {}", entrez_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeEntrezResource>(entrez_ident, entrez_file_name);

      auto const [iter, result] = resource_map.try_emplace(entrez_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Entrez Gene ident: {} to map (duplicate)", entrez_ident);

      }

    } else if (tree_type == PF7_SAMPLE_DATABASE_) {

      std::string Pf7Sample_ident;
      if (not sub_tree.getProperty(PF7_SAMPLE_IDENT_, Pf7Sample_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample Identifier.");
        continue;

      }

      std::string Pf7Sample_file_name;
      if (not sub_tree.getFileProperty(PF7_SAMPLE_FILE_, workDirectory(), Pf7Sample_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pf7 sample data file name information, ident: {}", Pf7Sample_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimePf7Resource>(Pf7Sample_ident, Pf7Sample_file_name);

      auto const [iter, result] = resource_map.try_emplace(Pf7Sample_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Entrez Gene ident: {} to map (duplicate)", Pf7Sample_ident);

      }

    } else if (tree_type == PUBMED_LIT_API_) {

      std::string pubmed_api_ident;
      if (not sub_tree.getProperty(PUBMED_LIT_IDENT_, pubmed_api_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Pubmed API Identifier.");
        continue;

      }

      std::string publication_cache_file;  // Note that the cache file(s) may not exist.
      if (not sub_tree.getFileCreateProperty(PUBMED_PUBLICATION_CACHE_, workDirectory(), publication_cache_file)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; Cannot create Pubmed publication API Cache file, ident: {}", pubmed_api_ident);
        continue;

      }

      std::string citation_cache_file;  // Note that the cache file(s) may not exist.
      if (not sub_tree.getFileCreateProperty(PUBMED_CITATION_CACHE_, workDirectory(), citation_cache_file)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; Cannot create Pubmed citation API Cache file, ident: {}", pubmed_api_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimePubmedAPIResource>( pubmed_api_ident,
                                                                                                              publication_cache_file,
                                                                                                              citation_cache_file);

      auto const [iter, result] = resource_map.try_emplace(pubmed_api_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Pubmed APi ident: {} to map (duplicate)", pubmed_api_ident);

      }

    } else if (tree_type == GENE_ID_DATABASE_) {

      std::string nomenclature_ident;
      if (not sub_tree.getProperty(GENE_ID_IDENT_, nomenclature_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No gene Nomenclature Identifier.");
        continue;

      }

      std::string nomenclature_file_name;
      if (not sub_tree.getFileProperty(GENE_ID_FILE_, workDirectory(), nomenclature_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Gene Nomenclature file name information, ident: {}", nomenclature_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeNomenclatureResource>(nomenclature_ident,
                                                                                                                nomenclature_file_name);

      auto const [iter, result] = resource_map.try_emplace(nomenclature_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Gene Nomenclature ident: {} to map (duplicate)", nomenclature_ident);

      }

    } else if (tree_type == AUX_ID_DATABASE_) {

      std::string genome_aux_ident;
      if (not sub_tree.getProperty(AUX_ID_IDENT_, genome_aux_ident)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Aux Identifier.");
        continue;

      }

      std::string genome_aux_file_name;
      if (not sub_tree.getFileProperty(AUX_ID_FILE_, workDirectory(), genome_aux_file_name)) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources; No Genome Aux file name information, ident: {}", genome_aux_ident);
        continue;

      }

      std::shared_ptr<const RuntimeResource> resource_ptr = std::make_shared<const RuntimeGenomeAuxResource>(genome_aux_ident, genome_aux_file_name);

      auto const [iter, result] = resource_map.try_emplace(genome_aux_ident, resource_ptr);
      if (not result) {

        ExecEnv::log().error("RuntimeProperties::getRuntimeResources, Could not add Genome Aux ident: {} to map (duplicate)", genome_aux_ident);

      }

    } else if (tree_type != COMMENT_){

      ExecEnv::log().warn("RuntimeProperties::getRuntimeResources, Unexpected resource: '{}'; ignored", tree_type);

    }

  } // For all resources.

  return resource_map;

}

