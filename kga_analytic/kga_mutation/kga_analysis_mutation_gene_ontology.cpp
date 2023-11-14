//
// Created by kellerberrin on 14/4/21.
//

#include "kel_exec_env.h"
#include "kga_analysis_mutation_gene_ontology.h"
#include "kgl_gene_cache_ontology.h"


namespace kga = kellerberrin::genome::analysis;
namespace kgl = kellerberrin::genome;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Define the PIMPL implementation to be the kgl::OntologyGeneCache object.
class kga::OntologyCache::OntologyGeneCache : public kgl::OntologyGeneCache {

public:
  using kgl::OntologyGeneCache::OntologyGeneCache;

};


kga::OntologyCache::OntologyCache( const std::vector<std::string>& gene_vector,
                                   const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                                   const std::shared_ptr<const kol::GoGraph>& go_graph_ptr) {

  initializeOntology(gene_vector, annotation_ptr, go_graph_ptr);

}


kga::OntologyCache::~OntologyCache() {


  gene_cache_ptr_ = nullptr;

}


void kga::OntologyCache::initializeOntology( const std::vector<std::string>& gene_vector,
                                             const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                                             const std::shared_ptr<const kol::GoGraph>& go_graph_ptr) {

  gene_cache_ptr_ = std::make_shared<const OntologyGeneCache>(gene_vector, annotation_ptr, go_graph_ptr);

}


double kga::OntologyCache::setSimilarityBP(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityBP(gene);

}


double kga::OntologyCache::setSimilarityMF(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityMF(gene);

}


double kga::OntologyCache::setSimilarityCC(const std::string& gene) const {

  return gene_cache_ptr_->setSimilarityCC(gene);

}

bool kga::OntologyCache::isTargetGene(const std::string& gene) const {

  return gene_cache_ptr_->isTargetGene(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxMFSim(const std::string& gene) const {

  return gene_cache_ptr_->maxMFSim(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxBPSim(const std::string& gene) const {

  return gene_cache_ptr_->maxBPSim(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxCCSim(const std::string& gene) const {

  return gene_cache_ptr_->maxCCSim(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxFunMFBPSim(const std::string& gene) const {

  return gene_cache_ptr_->maxFunMFBPSim(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxMFInfo(const std::string& gene) const {

  return gene_cache_ptr_->maxMFInfo(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxBPInfo(const std::string& gene) const {

  return gene_cache_ptr_->maxBPInfo(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxCCInfo(const std::string& gene) const {

  return gene_cache_ptr_->maxCCInfo(gene);

}

std::pair<std::string, double> kga::OntologyCache::maxFunMFBPInfo(const std::string& gene) const {

  return gene_cache_ptr_->maxFunMFBPInfo(gene);

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void kga::OntologyStats::processOntologyStats(const std::string& gene_id,
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

  // Similarity
  auto [max_MF_gene_sim, max_MF_sim] = ontology_cache.maxMFSim(gene_id);
  max_MF_sim_ = max_MF_sim;
  max_MF_gene_sim_ = max_MF_gene_sim;

  auto [max_BP_gene_sim, max_BP_sim] = ontology_cache.maxBPSim(gene_id);
  max_BP_sim_ = max_BP_sim;
  max_BP_gene_sim_ = max_BP_gene_sim;

  auto [max_CC_gene_sim, max_CC_sim] = ontology_cache.maxCCSim(gene_id);
  max_CC_sim_ = max_CC_sim;
  max_CC_gene_sim_ = max_CC_gene_sim;

  auto [max_FunMFBP_gene_sim, max_FunMFBP_sim] = ontology_cache.maxFunMFBPSim(gene_id);
  max_FunMFBP_sim_ = max_FunMFBP_sim;
  max_FunMFBP_gene_sim_ = max_FunMFBP_gene_sim;

  // Information
  auto [max_MF_gene_info, max_MF_info] = ontology_cache.maxMFInfo(gene_id);
  max_MF_info_ = max_MF_info;
  max_MF_gene_info_ = max_MF_gene_info;

  auto [max_BP_gene_info, max_BP_info] = ontology_cache.maxBPInfo(gene_id);
  max_BP_info_ = max_BP_info;
  max_BP_gene_info_ = max_BP_gene_info;

  auto [max_CC_gene_info, max_CC_info] = ontology_cache.maxCCInfo(gene_id);
  max_CC_info_ = max_CC_info;
  max_CC_gene_info_ = max_CC_gene_info;

  auto [max_FunMFBP_gene_info, max_FunMFBP_info] = ontology_cache.maxFunMFBPInfo(gene_id);
  max_FunMFBP_info_ = max_FunMFBP_info;
  max_FunMFBP_gene_info_ = max_FunMFBP_gene_info;




}




void kga::OntologyStats::writeOntology(std::ostream& out_file, char output_delimiter) const {

  out_file << score_BP_ << output_delimiter
      << score_MF_ << output_delimiter
      << score_CC_ << output_delimiter
      << max_score_ << output_delimiter
      << av_score_ << output_delimiter
      << (targetGene_ ? "Malaria" : "") << output_delimiter
      << max_MF_sim_ << output_delimiter
      << max_MF_gene_sim_ << output_delimiter
      << max_BP_sim_ << output_delimiter
      << max_BP_gene_sim_ << output_delimiter
      << max_CC_sim_ << output_delimiter
      << max_CC_gene_sim_ << output_delimiter
      << max_FunMFBP_sim_ << output_delimiter
      << max_FunMFBP_gene_sim_ << output_delimiter
      << max_MF_info_ << output_delimiter
      << max_MF_gene_info_ << output_delimiter
      << max_BP_info_ << output_delimiter
      << max_BP_gene_info_ << output_delimiter
      << max_CC_info_ << output_delimiter
      << max_CC_gene_info_ << output_delimiter
      << max_FunMFBP_info_ << output_delimiter
      << max_FunMFBP_gene_info_;

}

void kga::OntologyStats::writeOntologyHeader(std::ostream& out_file, char output_delimiter) const {

  out_file << "ScoreBP" << output_delimiter
           << "ScoreMF" << output_delimiter
           << "ScoreCC" << output_delimiter
           << "MaxScore" << output_delimiter
           << "AvScore" << output_delimiter
           << "TargetGene" << output_delimiter
           << "max_MF_sim" << output_delimiter
           << "max_MF_gene_sim" << output_delimiter
           << "max_BP_sim" << output_delimiter
           << "max_BP_gene_sim" << output_delimiter
           << "max_CC_sim" << output_delimiter
           << "max_CC_gene_sim" << output_delimiter
           << "max_FunMFBP_sim" << output_delimiter
           << "max_FunMFBP_gene_sim" << output_delimiter
           << "max_MF_info" << output_delimiter
           << "max_MF_gene_info" << output_delimiter
           << "max_BP_info" << output_delimiter
           << "max_BP_gene_info" << output_delimiter
           << "max_CC_info" << output_delimiter
           << "max_CC_gene_info" << output_delimiter
           << "max_FunMFBP_info" << output_delimiter
           << "max_FunMFBP_gene_info";


}
