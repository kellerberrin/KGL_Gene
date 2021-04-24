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

kgl::OntologyCache::OntologyCache(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  initializeOntology(ontology_db_ptr);

}


kgl::OntologyCache::~OntologyCache() {


  gene_cache_ptr_ = nullptr;

}


void kgl::OntologyCache::initializeOntology(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) {

  // Create a vector of Malaria genes from the map.
  std::vector<std::string> gene_vector;
  for (auto const&[go_gene, gene_name] : malaria_gene_map_) {

    gene_vector.push_back(go_gene);

  }


  gene_cache_ptr_ = std::make_shared<const OntologyGeneCache>(gene_vector, ontology_db_ptr);

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
