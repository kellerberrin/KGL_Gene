//
// Created by kellerberrin on 24/4/21.
//

#ifndef KOL_GENE_CACHE_ONTOLOGY_H
#define KOL_GENE_CACHE_ONTOLOGY_H

#include <fstream>
#include "kgl_ontology_database.h"
#include "kol_GoEnums.h"
#include "kol_AsymmetricSimilarityCache.h"
#include "kol_SimilarityLin.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityJiangConrath.h"
#include "kol_BestMatchAverageSetSimilarity.h"


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
                     const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr)
      : ontology_db_ptr_(ontology_db_ptr) {

    gene_set_ = kol::SetUtilities::convertVector(gene_vector);
    initializeOntology(ontology_db_ptr_);

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
  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_BP_;
  std::vector<std::string> all_go_terms_BP_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_BP_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_BP_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_MF_;
  std::vector<std::string> all_go_terms_MF_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_MF_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_MF_ptr_;

  kol::OntologySetType<std::string> gene_set_go_terms_CC_;
  std::vector<std::string> all_go_terms_CC_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_CC_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_CC_ptr_;


  void initializeOntology(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  kol::OntologySetType<std::string> getGeneSetGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                                       kol::GO::Ontology ontology);
  std::vector<std::string> getAllGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                                          kol::GO::Ontology ontology);
  std::shared_ptr<const kol::LinSimilarity> getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  std::shared_ptr<const kol::ResnikSimilarity> getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  std::shared_ptr<const kol::JiangConrathSimilarity> getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  // Cached GO Term and Gene Data.

  [[nodiscard]] std::vector<std::string> geneSetGOTermsBP() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_BP_.begin(), gene_set_go_terms_BP_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsBP() const { return all_go_terms_BP_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache& cacheBP() const { return *cache_BP_ptr_; }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsMF() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_MF_.begin(), gene_set_go_terms_MF_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsMF() const { return all_go_terms_MF_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache& cacheMF() const { return *cache_MF_ptr_; }

  [[nodiscard]] std::vector<std::string> geneSetGOTermsCC() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(gene_set_go_terms_CC_.begin(), gene_set_go_terms_CC_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsCC() const { return all_go_terms_CC_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache&  cacheCC() const { return *cache_CC_ptr_; }

};



} // namespace

#endif //KOL_GENE_ONTOLOGY_H
