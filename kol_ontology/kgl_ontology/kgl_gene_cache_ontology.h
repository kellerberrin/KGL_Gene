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
#include "kol_SimilarityLin.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityJiangConrath.h"
#include "kol_SetSimilarityBestMatchAverage.h"


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

private:

  kol::OntologySetType<std::string> gene_set_;
  std::shared_ptr<const kol::TermAnnotation> annotation_ptr_;
  std::shared_ptr<const kol::GoGraph> go_graph_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_BP_;
  std::vector<std::string> all_go_terms_BP_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_BP_ptr_;
  std::shared_ptr<const kol::SetSimilarityBestMatchAverage> set_sim_BP_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_MF_;
  std::vector<std::string> all_go_terms_MF_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_MF_ptr_;
  std::shared_ptr<const kol::SetSimilarityBestMatchAverage> set_sim_MF_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_CC_;
  std::vector<std::string> all_go_terms_CC_;
  std::shared_ptr<const kol::SimilarityCacheAsymmetric> cache_CC_ptr_;
  std::shared_ptr<const kol::SetSimilarityBestMatchAverage> set_sim_CC_ptr_;


  void initializeOntology();
  kol::OntologySetType<std::string> getGeneSetGOVector(kol::GO::Ontology ontology);
  std::vector<std::string> getAllGOVector(kol::GO::Ontology ontology);
  std::shared_ptr<const kol::SimilarityLin> getLinSimilarity();
  std::shared_ptr<const kol::SimilarityResnik> getResnikSimilarity();
  std::shared_ptr<const kol::SimilarityJIangConrath> getJiangSimilarity();
  // Cached GO Term and Gene Data.

  [[nodiscard]] std::vector<std::string> geneSetGOTermsBP() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_BP_.begin(), gene_set_go_terms_BP_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsBP() const { return all_go_terms_BP_; }
  [[nodiscard]] const kol::SimilarityCacheAsymmetric& cacheBP() const { return *cache_BP_ptr_; }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsMF() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_MF_.begin(), gene_set_go_terms_MF_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsMF() const { return all_go_terms_MF_; }
  [[nodiscard]] const kol::SimilarityCacheAsymmetric& cacheMF() const { return *cache_MF_ptr_; }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsCC() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_CC_.begin(), gene_set_go_terms_CC_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsCC() const { return all_go_terms_CC_; }
  [[nodiscard]] const kol::SimilarityCacheAsymmetric&  cacheCC() const { return *cache_CC_ptr_; }

};



} // namespace

#endif //KOL_GENE_ONTOLOGY_H
