//
// Created by kellerberrin on 24/4/21.
//

#ifndef KOL_GENE_CACHE_ONTOLOGY_H
#define KOL_GENE_CACHE_ONTOLOGY_H

#include <fstream>
#include "kol_SetUtilities.h"
#include "kgl_ontology_database.h"
#include "kol_GoEnums.h"
#include "kol_SimilarityCacheAsymmetric.h"
#include "kol_SimilarityInterface.h"
#include "kol_SetSimilarityInterface.h"
#include "kol_OntologyFactory.h"


namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class OntologyGeneCache {

public:

  OntologyGeneCache( const std::vector<std::string>& gene_vector,
                     const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                     const std::shared_ptr<const kol::GoGraph>& go_graph_ptr)
      : annotation_ptr_(annotation_ptr), go_graph_ptr_(go_graph_ptr) {

    gene_set_ = kol::SetUtilities::convertVector(gene_vector);
    initializeOntology();

  }
  OntologyGeneCache(const OntologyGeneCache &) = delete;
  ~OntologyGeneCache() = default;

  OntologyGeneCache &operator=(const OntologyGeneCache &) = delete;


  [[nodiscard]] double setSimilarityBP(const std::string& gene) const;
  [[nodiscard]] double setSimilarityMF(const std::string& gene) const;
  [[nodiscard]] double setSimilarityCC(const std::string& gene) const;
  [[nodiscard]] bool isTargetGene(const std::string& gene) const { return gene_set_.contains(gene); }
  [[nodiscard]] std::pair<std::string, double> maxMF(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxBP(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxCC(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxFunSim(const std::string& gene) const;

private:

  // Options are:  ANCESTOR, CONTENT, CONTENT_DAG, COUTOGRASM, COUTOGRASMADJUSTED, EXCLUSIVEINHERITED, FRONTIER
  static const constexpr kol::InformationContentType information_type_ = kol::InformationContentType::CONTENT;
  // Options are: RESNIK, LIN, RELEVANCE, JIANG, PEKARSTAAB
  static const constexpr kol::SimilarityType similarity_type_ = kol::SimilarityType::LIN;
  // Options are: GENTLEMANSIMUI, JACCARD.  Stand-alone.
  //  MAZANDUSIMDIC, MAZANDUSIMUIC, PESQUITASIMGIC.  Information source.
  //  BESTMATCHAVERAGE, AVERAGEBESTMATCH, ALLPAIRSMAX, ALLPAIRSAVERAGE.  Similarity measure.
  static const constexpr kol::SetSimilarityType setsimilarity_type_ = kol::SetSimilarityType::AVERAGEBESTMATCH;

  kol::OntologySetType<std::string> gene_set_;
  std::shared_ptr<const kol::TermAnnotation> annotation_ptr_;
  std::shared_ptr<const kol::GoGraph> go_graph_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_BP_;
  std::vector<std::string> all_go_terms_BP_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_BP_ptr_;
  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_BP_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_MF_;
  std::vector<std::string> all_go_terms_MF_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_MF_ptr_;
  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_MF_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_CC_;
  std::vector<std::string> all_go_terms_CC_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_CC_ptr_;
  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_CC_ptr_;


  void initializeOntology();
  kol::OntologySetType<std::string> getGeneSetGOVector(kol::GO::Ontology ontology);
  std::vector<std::string> getAllGOVector(kol::GO::Ontology ontology);
  std::shared_ptr<const kol::SimilarityInterface> getSimilarity();
  // Cached GO Term and Gene Data.

  [[nodiscard]] std::vector<std::string> geneSetGOTermsBP() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_BP_.begin(), gene_set_go_terms_BP_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsMF() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_MF_.begin(), gene_set_go_terms_MF_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsCC() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_CC_.begin(), gene_set_go_terms_CC_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }

};



} // namespace

#endif //KOL_GENE_ONTOLOGY_H
