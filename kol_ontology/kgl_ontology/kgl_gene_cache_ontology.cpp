//
// Created by kellerberrin on 24/4/21.
//

#include "kgl_gene_cache_ontology.h"
#include "kgl_ontology_database_test.h"
#include "kol_InformationContent.h"
#include "kel_exec_env.h"

namespace kol = kellerberrin::ontology;
namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::OntologyGeneCache::initializeOntology(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  // Remove when testing complete.
  OntologyDatabaseTest database_test(ontology_db_ptr);
  database_test.performTests();

  auto sim_calc_ptr = getLinSimilarity(ontology_db_ptr);

  // Biological process
  gene_set_go_terms_BP_ = getGeneSetGOVector(ontology_db_ptr, kol::GO::Ontology::BIOLOGICAL_PROCESS);
  all_go_terms_BP_ = getAllGOVector(ontology_db_ptr, kol::GO::Ontology::BIOLOGICAL_PROCESS);
  cache_BP_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsBP(),
                                                                         all_go_terms_BP_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::BIOLOGICAL_PROCESS);
  set_sim_BP_ptr_ = std::make_shared<const kol::SetSimilarityBestMatchAverage>(cache_BP_ptr_);
  ExecEnv::log().info("BP GO cache created, All BP terms (columns): {}, Gene set BP terms (rows): {}", cache_BP_ptr_->columns(), cache_BP_ptr_->rows());

  // Molecular Function
  gene_set_go_terms_MF_ = getGeneSetGOVector(ontology_db_ptr, kol::GO::Ontology::MOLECULAR_FUNCTION);
  all_go_terms_MF_ = getAllGOVector(ontology_db_ptr, kol::GO::Ontology::MOLECULAR_FUNCTION);
  cache_MF_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsMF(),
                                                                         all_go_terms_MF_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::MOLECULAR_FUNCTION);
  set_sim_MF_ptr_ = std::make_shared<const kol::SetSimilarityBestMatchAverage>(cache_MF_ptr_);
  ExecEnv::log().info("MF GO cache created, All MF terms (columns): {}, Gene set MF terms (rows): {}", cache_MF_ptr_->columns(), cache_MF_ptr_->rows());

  // Cellular Component
  gene_set_go_terms_CC_ = getGeneSetGOVector(ontology_db_ptr, kol::GO::Ontology::CELLULAR_COMPONENT);
  all_go_terms_CC_ =  getAllGOVector(ontology_db_ptr, kol::GO::Ontology::CELLULAR_COMPONENT);
  cache_CC_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsCC(),
                                                                         all_go_terms_CC_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::CELLULAR_COMPONENT);
  set_sim_CC_ptr_ = std::make_shared<const kol::SetSimilarityBestMatchAverage>(cache_CC_ptr_);
  ExecEnv::log().info("CC GO cache created, All CC terms (columns): {}, Gene set CC terms (rows): {}",  cache_CC_ptr_->columns(), cache_CC_ptr_->rows());

}

kol::OntologySetType<std::string> kgl::OntologyGeneCache::getGeneSetGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                                                             kol::GO::Ontology ontology) {

  const kol::TermAnnotation& annotation = *ontology_db_ptr->annotation();

  // Only unique GO terms.
  kol::OntologySetType<std::string> unique_GO_terms;
  for (auto const& gene : gene_set_) {

    for (auto const& go_term : annotation.getGoTermsForGeneByOntology(gene, ontology)) {

      unique_GO_terms.insert(go_term);

    }

  }

  return unique_GO_terms;

}

std::vector<std::string> kgl::OntologyGeneCache::getAllGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr, kol::GO::Ontology ontology) {

  const kol::TermAnnotation& annotation = *ontology_db_ptr->annotation();

  return annotation.getOntologyTerms(ontology);

}


std::shared_ptr<const kol::SimilarityLin> kgl::OntologyGeneCache::getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {


  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityLin>(ic_map_ptr);

}


std::shared_ptr<const kol::SimilarityResnik> kgl::OntologyGeneCache::getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityResnik>(ic_map_ptr);

}

std::shared_ptr<const kol::SimilarityJIangConrath> kgl::OntologyGeneCache::getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityJIangConrath>(ic_map_ptr);

}


double kgl::OntologyGeneCache::setSimilarityBP(const std::string& gene) const {

  const kol::TermAnnotation& annotation = *ontology_db_ptr_->annotation();

  auto gene_BP_terms = annotation.getGoTermsForGeneBP(gene);

  return set_sim_BP_ptr_->calculateSimilarity(gene_set_go_terms_BP_, gene_BP_terms);

}

double kgl::OntologyGeneCache::setSimilarityMF(const std::string& gene) const {

  const kol::TermAnnotation& annotation = *ontology_db_ptr_->annotation();

  auto gene_MF_terms = annotation.getGoTermsForGeneMF(gene);

  return set_sim_MF_ptr_->calculateSimilarity(gene_set_go_terms_MF_, gene_MF_terms);

}

double kgl::OntologyGeneCache::setSimilarityCC(const std::string& gene) const {

  const kol::TermAnnotation& annotation = *ontology_db_ptr_->annotation();

  auto gene_CC_terms = annotation.getGoTermsForGeneCC(gene);

  return set_sim_CC_ptr_->calculateSimilarity(gene_set_go_terms_CC_, gene_CC_terms);

}



