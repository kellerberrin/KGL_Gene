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
    initializeSimilarity();
    initializeInformation();

  }
  OntologyGeneCache(const OntologyGeneCache &) = delete;
  ~OntologyGeneCache() = default;

  OntologyGeneCache &operator=(const OntologyGeneCache &) = delete;


  [[nodiscard]] double setSimilarityBP(const std::string& gene) const;
  [[nodiscard]] double setSimilarityMF(const std::string& gene) const;
  [[nodiscard]] double setSimilarityCC(const std::string& gene) const;

  [[nodiscard]] bool isTargetGene(const std::string& gene) const { return gene_set_.contains(gene); }

  [[nodiscard]] std::pair<std::string, double> maxMFSim(const std::string& check_gene) const { return maxMF(check_gene, set_sim_MF_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxBPSim(const std::string& check_gene) const { return maxBP(check_gene, set_sim_BP_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxCCSim(const std::string& check_gene) const { return maxCC(check_gene, set_sim_CC_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxFunMFBPSim(const std::string& check_gene) const { return maxFunMFBP(check_gene, set_sim_MF_ptr_, set_sim_BP_ptr_); }

  [[nodiscard]] std::pair<std::string, double> maxMFInfo(const std::string& check_gene) const { return maxMF(check_gene, set_info_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxBPInfo(const std::string& check_gene) const { return maxBP(check_gene, set_info_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxCCInfo(const std::string& check_gene) const { return maxCC(check_gene, set_info_ptr_); }
  [[nodiscard]] std::pair<std::string, double> maxFunMFBPInfo(const std::string& check_gene) const { return maxFunMFBP(check_gene, set_info_ptr_, set_info_ptr_); }

private:

  // Input objects.
  kol::OntologySetType<std::string> gene_set_;
  std::shared_ptr<const kol::TermAnnotation> annotation_ptr_;
  std::shared_ptr<const kol::GoGraph> go_graph_ptr_;

  // GO term sets.
  kol::OntologySetType<std::string> gene_set_go_terms_BP_;
  std::vector<std::string> all_go_terms_BP_;

  kol::OntologySetType<std::string> gene_set_go_terms_MF_;
  std::vector<std::string> all_go_terms_MF_;

  kol::OntologySetType<std::string> gene_set_go_terms_CC_;
  std::vector<std::string> all_go_terms_CC_;

  // Similarity based members
  // Options are:  ANCESTOR, CONTENT, CONTENT_DAG, COUTOGRASM, COUTOGRASMADJUSTED, EXCLUSIVEINHERITED, FRONTIER
  static const constexpr kol::InformationContentType similarity_info_type_ = kol::InformationContentType::CONTENT;
  // Options are: RESNIK, LIN, RELEVANCE, JIANG, PEKARSTAAB
  static const constexpr kol::SimilarityType similarity_type_ = kol::SimilarityType::LIN;
  //  Options are: BESTMATCHAVERAGE, AVERAGEBESTMATCH, ALLPAIRSMAX, ALLPAIRSAVERAGE.  Similarity measure.
  static const constexpr kol::SetSimilarityType setsimilarity_type_ = kol::SetSimilarityType::BESTMATCHAVERAGE;

  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_BP_ptr_;
  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_MF_ptr_;
  std::shared_ptr<const kol::SetSimilarityInterface> set_sim_CC_ptr_;

  // Information based members.
  // Options are:  ANCESTOR, CONTENT, CONTENT_DAG, COUTOGRASM, COUTOGRASMADJUSTED, EXCLUSIVEINHERITED, FRONTIER
  static const constexpr kol::InformationContentType information_type_ = kol::InformationContentType::CONTENT;
  //  Options are: MAZANDUSIMDIC, MAZANDUSIMUIC, PESQUITASIMGIC.  Information source.
  static const constexpr kol::SetSimilarityType setinfo_type_ = kol::SetSimilarityType::MAZANDUSIMUIC;

  std::shared_ptr<const kol::SetSimilarityInterface> set_info_ptr_;

  void initializeSimilarity();
  void initializeInformation();

  [[nodiscard]] kol::OntologySetType<std::string> getGeneSetGOVector(kol::GO::Ontology ontology) const;
  [[nodiscard]] std::vector<std::string> getAllGOVector(kol::GO::Ontology ontology) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getSimilarity() const;
  // Cached GO Term and Gene Data.
  [[nodiscard]] std::vector<std::string> geneSetGOTermsBP() const;
  [[nodiscard]] std::vector<std::string> geneSetGOTermsMF() const;
  [[nodiscard]] std::vector<std::string> geneSetGOTermsCC() const;
  // Set Similarity Dependent calculations.
  [[nodiscard]] std::pair<std::string, double> maxMF(const std::string& check_gene,
                                                     const std::shared_ptr<const kol::SetSimilarityInterface>& set_sim_MF_ptr) const;
  [[nodiscard]] std::pair<std::string, double> maxBP(const std::string& check_gene,
                                                     const std::shared_ptr<const kol::SetSimilarityInterface>& set_sim_BP_ptr) const;
  [[nodiscard]] std::pair<std::string, double> maxCC(const std::string& check_gene,
                                                     const std::shared_ptr<const kol::SetSimilarityInterface>& set_sim_CC_ptr) const;
  [[nodiscard]] std::pair<std::string, double> maxFunMFBP(const std::string& check_gene,
                                                          const std::shared_ptr<const kol::SetSimilarityInterface>& set_sim_MF_ptr,
                                                          const std::shared_ptr<const kol::SetSimilarityInterface>& set_sim_BP_ptr) const;

};



} // namespace

#endif //KOL_GENE_ONTOLOGY_H
