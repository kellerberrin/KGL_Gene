//
// Created by kellerberrin on 14/4/21.
//

#include "kel_exec_env.h"
#include "kgl_analysis_mutation_gene_ontology.h"


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::OntologyCache::initializeOntology(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  // Create a vector of Malaria genes from the map.
  for (auto const& [Go_gene, Gene_name] : malaria_gene_map_) {

    malaria_gene_vector_.push_back(Go_gene);

  }

  auto sim_calc_ptr = getResnikSimilarity(ontology_db_ptr);

  // Biological process
  malaria_go_terms_BP_ = getMalariaGOVector(ontology_db_ptr, kol::GO::Ontology::BIOLOGICAL_PROCESS);
  all_go_terms_BP_ = getAllGOVector(ontology_db_ptr, kol::GO::Ontology::BIOLOGICAL_PROCESS);
  cache_BP_ptr_ = std::make_shared<const kol::AsymmetricSimilarityCache>(malariaGOTermsBP(),
                                                                         all_go_terms_BP_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::BIOLOGICAL_PROCESS);
  set_sim_BP_ptr_ = std::make_shared<const kol::BestMatchAverageSetSimilarity>(cache_BP_ptr_);

  // Molecular Function
  malaria_go_terms_MF_ = getMalariaGOVector(ontology_db_ptr, kol::GO::Ontology::MOLECULAR_FUNCTION);
  all_go_terms_MF_ = getAllGOVector(ontology_db_ptr, kol::GO::Ontology::MOLECULAR_FUNCTION);
  cache_MF_ptr_ = std::make_shared<const kol::AsymmetricSimilarityCache>(malariaGOTermsMF(),
                                                                         all_go_terms_MF_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::MOLECULAR_FUNCTION);
  set_sim_MF_ptr_ = std::make_shared<const kol::BestMatchAverageSetSimilarity>(cache_MF_ptr_);

  // Cellular Component
  malaria_go_terms_CC_ = getMalariaGOVector(ontology_db_ptr, kol::GO::Ontology::CELLULAR_COMPONENT);
  all_go_terms_CC_ =  getAllGOVector(ontology_db_ptr, kol::GO::Ontology::CELLULAR_COMPONENT);
  cache_CC_ptr_ = std::make_shared<const kol::AsymmetricSimilarityCache>(malariaGOTermsCC(),
                                                                         all_go_terms_CC_,
                                                                         sim_calc_ptr,
                                                                         kol::GO::Ontology::CELLULAR_COMPONENT);
  set_sim_CC_ptr_ = std::make_shared<const kol::BestMatchAverageSetSimilarity>(cache_CC_ptr_);

}

kol::OntologySetType<std::string> kgl::OntologyCache::getMalariaGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr, kol::GO::Ontology ontology) {

  const kol::AnnotationData& annotation = *ontology_db_ptr->annotation();
  const kol::GoGraph& go_graph = *ontology_db_ptr->goGraph();

  // Only unique GO terms.
  kol::OntologySetType<std::string> unique_GO_terms;
  for (auto const& gene : malaria_gene_vector_) {

    for (auto const& go_term : annotation.getGoTermsForGeneByOntology(gene, ontology, go_graph)) {

      unique_GO_terms.insert(go_term);

    }

  }

  return unique_GO_terms;

}

std::vector<std::string> kgl::OntologyCache::getAllGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr, kol::GO::Ontology ontology) {

  const kol::AnnotationData& annotation = *ontology_db_ptr->annotation();
  const kol::GoGraph& go_graph = *ontology_db_ptr->goGraph();

  return annotation.getOntologyTerms(go_graph, ontology);

}


std::shared_ptr<const kol::LinSimilarity> kgl::OntologyCache::getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  std::shared_ptr<const kol::TermInformationContentMap> info_map_ptr(std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr->goGraph(),
                                                                                                                            ontology_db_ptr->annotation()));
  return std::make_shared<const kol::LinSimilarity>(ontology_db_ptr->goGraph(), info_map_ptr);

}


std::shared_ptr<const kol::ResnikSimilarity> kgl::OntologyCache::getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  std::shared_ptr<const kol::TermInformationContentMap> info_map_ptr(std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr->goGraph(),
                                                                                                                            ontology_db_ptr->annotation()));
  return std::make_shared<const kol::ResnikSimilarity>(ontology_db_ptr->goGraph(), info_map_ptr);

}


[[nodiscard]] double kgl::OntologyCache::setSimilarityBP(const std::string& gene) const {

  const kol::AnnotationData& annotation = *ontology_db_ptr_->annotation();
  const kol::GoGraph& go_graph = *ontology_db_ptr_->goGraph();

  auto gene_BP_terms = annotation.getGoTermsForGeneBP(gene, go_graph);

  return set_sim_BP_ptr_->calculateSimilarity(malaria_go_terms_BP_, gene_BP_terms);

}

[[nodiscard]] double kgl::OntologyCache::setSimilarityMF(const std::string& gene) const {

  const kol::AnnotationData& annotation = *ontology_db_ptr_->annotation();
  const kol::GoGraph& go_graph = *ontology_db_ptr_->goGraph();

  auto gene_MF_terms = annotation.getGoTermsForGeneMF(gene, go_graph);

  return set_sim_MF_ptr_->calculateSimilarity(malaria_go_terms_MF_, gene_MF_terms);

}

[[nodiscard]] double kgl::OntologyCache::setSimilarityCC(const std::string& gene) const {

  const kol::AnnotationData& annotation = *ontology_db_ptr_->annotation();
  const kol::GoGraph& go_graph = *ontology_db_ptr_->goGraph();

  auto gene_CC_terms = annotation.getGoTermsForGeneCC(gene, go_graph);

  return set_sim_CC_ptr_->calculateSimilarity(malaria_go_terms_CC_, gene_CC_terms);

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::OntologyStats::processOntologyStats(const GeneCharacteristic& gene_info,
                                              const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                              const OntologyCache& ontology_cache) {

  score_BP_ = ontology_cache.setSimilarityBP(gene_info.gafId());
  score_MF_ = ontology_cache.setSimilarityMF(gene_info.gafId());
  score_CC_ = ontology_cache.setSimilarityCC(gene_info.gafId());
  max_score_ = std::max(score_BP_, std::max(score_MF_, score_CC_));
  double sum {0.0}, count{0.0};
  if (score_BP_ > 0.0) {

    sum += score_BP_;
    count += 1.0;

  }
  if (score_MF_ > 0.0) {

    sum += score_MF_;
    count += 1.0;

  }
  if (score_CC_ > 0.0) {

    sum += score_CC_;
    count += 1.0;

  }

  if (count > 0.0) {

    av_score_ = sum / count;

  }

  targetGene_ = ontology_cache.isTargetGene(gene_info.gafId());

}



void kgl::OntologyStats::writeOntology(std::ostream& out_file, char output_delimiter) const {

  out_file << score_BP_ << output_delimiter
           << score_MF_ << output_delimiter
           << score_CC_ << output_delimiter
           << max_score_ << output_delimiter
           << av_score_ << output_delimiter
           << (targetGene_ ? "Malaria" : "");

}

void kgl::OntologyStats::writeOntologyHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "ScoreBP" << output_delimiter
           << "ScoreMF" << output_delimiter
           << "ScoreCC" << output_delimiter
           << "MaxScore" << output_delimiter
           << "AvScore" << output_delimiter
           << "TargetGene";

}
