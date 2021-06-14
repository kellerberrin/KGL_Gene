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



void kgl::OntologyGeneCache::initializeOntology() {


  auto sim_calc_ptr = getSimilarity();

  // Biological process
  gene_set_go_terms_BP_ = getGeneSetGOVector(kol::GO::Ontology::BIOLOGICAL_PROCESS);
  all_go_terms_BP_ = getAllGOVector(kol::GO::Ontology::BIOLOGICAL_PROCESS);
  cache_BP_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsBP(),
                                                                         all_go_terms_BP_,
                                                                         sim_calc_ptr);
  set_sim_BP_ptr_ = kol::OntologyFactory::createSetSimilarity(cache_BP_ptr_, setsimilarity_type_);
  ExecEnv::log().info("BP GO cache created, All BP terms (columns): {}, Gene set BP terms (rows): {}", cache_BP_ptr_->columns(), cache_BP_ptr_->rows());

  // Molecular Function
  gene_set_go_terms_MF_ = getGeneSetGOVector(kol::GO::Ontology::MOLECULAR_FUNCTION);
  all_go_terms_MF_ = getAllGOVector(kol::GO::Ontology::MOLECULAR_FUNCTION);
  cache_MF_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsMF(),
                                                                         all_go_terms_MF_,
                                                                         sim_calc_ptr);
  set_sim_MF_ptr_ = kol::OntologyFactory::createSetSimilarity(cache_MF_ptr_, setsimilarity_type_);
  ExecEnv::log().info("MF GO cache created, All MF terms (columns): {}, Gene set MF terms (rows): {}", cache_MF_ptr_->columns(), cache_MF_ptr_->rows());

  // Cellular Component
  gene_set_go_terms_CC_ = getGeneSetGOVector(kol::GO::Ontology::CELLULAR_COMPONENT);
  all_go_terms_CC_ =  getAllGOVector(kol::GO::Ontology::CELLULAR_COMPONENT);
  cache_CC_ptr_ = std::make_shared<const kol::SimilarityCacheAsymmetric>(geneSetGOTermsCC(),
                                                                         all_go_terms_CC_,
                                                                         sim_calc_ptr);
  set_sim_CC_ptr_ = kol::OntologyFactory::createSetSimilarity(cache_CC_ptr_, setsimilarity_type_);
  ExecEnv::log().info("CC GO cache created, All CC terms (columns): {}, Gene set CC terms (rows): {}",  cache_CC_ptr_->columns(), cache_CC_ptr_->rows());

}

kol::OntologySetType<std::string> kgl::OntologyGeneCache::getGeneSetGOVector(kol::GO::Ontology ontology) {

  // Only unique GO terms.
  kol::OntologySetType<std::string> unique_GO_terms;
  for (auto const& gene : gene_set_) {

    for (auto const& go_term : annotation_ptr_->getGoTermsForGeneByOntology(gene, ontology)) {

      unique_GO_terms.insert(go_term);

    }

  }

  return unique_GO_terms;

}

std::vector<std::string> kgl::OntologyGeneCache::getAllGOVector(kol::GO::Ontology ontology) {

  return annotation_ptr_->getOntologyTerms(ontology);

}


std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyGeneCache::getSimilarity() {

  return kol::OntologyFactory::createSimilarity( annotation_ptr_, go_graph_ptr_, similarity_type_, information_type_);

}


double kgl::OntologyGeneCache::setSimilarityBP(const std::string& gene) const {

  auto gene_BP_terms = annotation_ptr_->getGoTermsForGeneBP(gene);

  return set_sim_BP_ptr_->calculateSimilarity(gene_set_go_terms_BP_, gene_BP_terms);

}

double kgl::OntologyGeneCache::setSimilarityMF(const std::string& gene) const {

  auto gene_MF_terms = annotation_ptr_->getGoTermsForGeneMF(gene);

  return set_sim_MF_ptr_->calculateSimilarity(gene_set_go_terms_MF_, gene_MF_terms);

}

double kgl::OntologyGeneCache::setSimilarityCC(const std::string& gene) const {

  auto gene_CC_terms = annotation_ptr_->getGoTermsForGeneCC(gene);

  return set_sim_CC_ptr_->calculateSimilarity(gene_set_go_terms_CC_, gene_CC_terms);

}

std::pair<std::string, double> kgl::OntologyGeneCache::maxMF(const std::string& gene) const {

  auto gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::MOLECULAR_FUNCTION);

  double max_similarity{0.0};
  std::string max_ref_gene;
  for (auto const& ref_gene : gene_set_) {

    auto ref_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::MOLECULAR_FUNCTION);

    double similarity = set_sim_MF_ptr_->calculateSimilarity(gene_go_set, ref_gene_go_set);

    if (similarity > max_similarity) {

      max_similarity = similarity;
      max_ref_gene = ref_gene;

    }

  }

  return std::pair<std::string, double>(max_ref_gene, max_similarity);

}


std::pair<std::string, double> kgl::OntologyGeneCache::maxBP(const std::string& gene) const {

  auto gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::BIOLOGICAL_PROCESS);

  double max_similarity{0.0};
  std::string max_ref_gene;
  for (auto const& ref_gene : gene_set_) {

    auto ref_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::BIOLOGICAL_PROCESS);

    double similarity = set_sim_BP_ptr_->calculateSimilarity(ref_gene_go_set, gene_go_set);

    if (similarity > max_similarity) {

      max_similarity = similarity;
      max_ref_gene = ref_gene;

    }

  }

  return std::pair<std::string, double>(max_ref_gene, max_similarity);

}


std::pair<std::string, double> kgl::OntologyGeneCache::maxCC(const std::string& gene) const {

  auto gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::CELLULAR_COMPONENT);

  double max_similarity{0.0};
  std::string max_ref_gene;
  for (auto const& ref_gene : gene_set_) {

    auto ref_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::CELLULAR_COMPONENT);

    double similarity = set_sim_CC_ptr_->calculateSimilarity(ref_gene_go_set, gene_go_set);

    if (similarity > max_similarity) {

      max_similarity = similarity;
      max_ref_gene = ref_gene;

    }

  }

  return std::pair<std::string, double>(max_ref_gene, max_similarity);

}


std::pair<std::string, double> kgl::OntologyGeneCache::maxFunSim(const std::string& gene) const {


  auto mf_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::MOLECULAR_FUNCTION);
  auto bp_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::BIOLOGICAL_PROCESS);

  double max_similarity{0.0};
  std::string max_ref_gene;
  for (auto const& ref_gene : gene_set_) {

    auto mf_ref_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::MOLECULAR_FUNCTION);
    auto bp_ref_gene_go_set = annotation_ptr_->getGoTermsForGeneByOntology(gene, kol::GO::Ontology::BIOLOGICAL_PROCESS);

    double mf_similarity = set_sim_MF_ptr_->calculateSimilarity(mf_ref_gene_go_set, mf_gene_go_set);
    double bp_similarity = set_sim_BP_ptr_->calculateSimilarity(bp_ref_gene_go_set, bp_gene_go_set);

    double similarity = ((mf_similarity * mf_similarity) + (bp_similarity * bp_similarity)) / 2.0;

    if (similarity > max_similarity) {

      max_similarity = similarity;
      max_ref_gene = ref_gene;

    }

  }

  return std::pair<std::string, double>(max_ref_gene, max_similarity);

}


