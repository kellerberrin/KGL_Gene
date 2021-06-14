//
// Created by kellerberrin on 14/4/21.
//

#include "kel_exec_env.h"
#include "kgl_analysis_mutation_gene_ontology.h"
#include "kgl_gene_cache_ontology.h"


namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::OntologyCache::OntologyCache( const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                                   const std::shared_ptr<const kol::GoGraph>& go_graph_ptr) {

  initializeOntology(annotation_ptr, go_graph_ptr);

}


kgl::OntologyCache::~OntologyCache() {


  gene_cache_ptr_ = nullptr;

}


void kgl::OntologyCache::initializeOntology( const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                                             const std::shared_ptr<const kol::GoGraph>& go_graph_ptr) {

  // Create a vector of Malaria genes from the map.
  std::vector<std::string> gene_vector;
  for (auto const&[go_gene, gene_name] : malaria_gene_map_) {

    gene_vector.push_back(gene_name);

  }

  gene_cache_ptr_ = std::make_shared<const OntologyGeneCache>(gene_vector, annotation_ptr, go_graph_ptr);

}


double kgl::OntologyCache::setSimilarityBP(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityBP(gene);

}


double kgl::OntologyCache::setSimilarityMF(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityMF(gene);

}


double kgl::OntologyCache::setSimilarityCC(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityCC(gene);

}

bool kgl::OntologyCache::isTargetGene(const std::string& gene) const {

  return gene_cache_ptr_->isTargetGene(gene);

}

std::pair<std::string, double> kgl::OntologyCache::maxMF(const std::string& gene) const {

  return gene_cache_ptr_->maxMF(gene);

}

std::pair<std::string, double> kgl::OntologyCache::maxBP(const std::string& gene) const {

  return gene_cache_ptr_->maxBP(gene);

}

std::pair<std::string, double> kgl::OntologyCache::maxCC(const std::string& gene) const {

  return gene_cache_ptr_->maxCC(gene);

}

std::pair<std::string, double> kgl::OntologyCache::maxFunSim(const std::string& gene) const {

  return gene_cache_ptr_->maxFunSim(gene);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kgl::OntologyStats::processOntologyStats(const std::string& gene_id,
                                              const OntologyCache& ontology_cache) {

  score_BP_ = ontology_cache.setSimilarityBP(gene_id);
  score_MF_ = ontology_cache.setSimilarityMF(gene_id);
  score_CC_ = ontology_cache.setSimilarityCC(gene_id);
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

  targetGene_ = ontology_cache.isTargetGene(gene_id);

  auto [max_MF_gene, max_MF] = ontology_cache.maxMF(gene_id);
  max_MF_ = max_MF;
  max_MF_gene_ = max_MF_gene;

//  auto [max_BP_gene, max_BP] = ontology_cache.maxBP(gene_id);
//  max_BP_ = max_BP;
//  max_BP_gene_ = max_BP_gene;

//  auto [max_CC_gene, max_CC] = ontology_cache.maxCC(gene_id);
//  max_CC_ = max_CC;
//  max_CC_gene_ = max_CC_gene;

//  auto [max_FunSim_gene, max_FunSim] = ontology_cache.maxFunSim(gene_id);
//  max_FunSim_ = max_FunSim;
//  max_FunSim_gene_ = max_FunSim_gene;

}




void kgl::OntologyStats::writeOntology(std::ostream& out_file, char output_delimiter) const {

  out_file << score_BP_ << output_delimiter
           << score_MF_ << output_delimiter
           << score_CC_ << output_delimiter
           << max_score_ << output_delimiter
           << av_score_ << output_delimiter
           << (targetGene_ ? "Malaria" : "") << output_delimiter
           << max_MF_ << output_delimiter
           << max_MF_gene_ << output_delimiter
           << max_BP_ << output_delimiter
           << max_BP_gene_ << output_delimiter
           << max_CC_ << output_delimiter
           << max_CC_gene_ << output_delimiter
           << max_FunSim_ << output_delimiter
           << max_FunSim_gene_;

}

void kgl::OntologyStats::writeOntologyHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "ScoreBP" << output_delimiter
           << "ScoreMF" << output_delimiter
           << "ScoreCC" << output_delimiter
           << "MaxScore" << output_delimiter
           << "AvScore" << output_delimiter
           << "TargetGene" << output_delimiter
           << "max_MF" << output_delimiter
           << "max_MF_gene" << output_delimiter
           << "max_BP" << output_delimiter
           << "max_BP_gene" << output_delimiter
           << "max_CC" << output_delimiter
           << "max_CC_gene" << output_delimiter
           << "max_FunSim" << output_delimiter
           << "max_FunSim_gene";


}
