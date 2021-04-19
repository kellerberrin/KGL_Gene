//
// Created by kellerberrin on 14/4/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
#define KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H

#include <fstream>
#include "kol_OntologyDatabase.h"
#include "kol_GoEnums.h"
#include "kol_AsymmetricSimilarityCache.h"
#include "kol_LinSimilarity.h"
#include "kol_ResnikSimilarity.h"
#include "kol_BestMatchAverageSetSimilarity.h"
#include "kgl_analysis_mutation_gene_stats.h"

namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class OntologyCache {

public:

  explicit OntologyCache(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr)
  : ontology_db_ptr_(ontology_db_ptr) { initializeOntology(ontology_db_ptr_); }
  OntologyCache(const OntologyCache &) = delete;
  ~OntologyCache() = default;

  OntologyCache &operator=(const OntologyCache &) = delete;


  [[nodiscard]] double setSimilarityBP(const std::string& gene) const;
  [[nodiscard]] double setSimilarityMF(const std::string& gene) const;
  [[nodiscard]] double setSimilarityCC(const std::string& gene) const;
  [[nodiscard]] bool isTargetGene(const std::string& gene) const { return malaria_gene_map_.contains(gene); }

private:

  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;
  std::vector<std::string> malaria_gene_vector_;

  kol::OntologySetType<std::string> malaria_go_terms_BP_;
  std::vector<std::string> all_go_terms_BP_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_BP_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_BP_ptr_;

  kol::OntologySetType<std::string> malaria_go_terms_MF_;
  std::vector<std::string> all_go_terms_MF_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_MF_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_MF_ptr_;

  kol::OntologySetType<std::string> malaria_go_terms_CC_;
  std::vector<std::string> all_go_terms_CC_;
  std::shared_ptr<const kol::AsymmetricSimilarityCache> cache_CC_ptr_;
  std::shared_ptr<const kol::BestMatchAverageSetSimilarity> set_sim_CC_ptr_;

  inline static const std::map<std::string, std::string> malaria_gene_map_ {
      { "P16671", "CD36"}, { "P06028","GYPB"}, { "P12318", "FCGR2A"}, { "P31994", "FCGR2B"}, { "P05362", "ICAM1"},
      {  "O14931", "NCR3"}, { "P68871", "HBB"}, { "P35228", "NOS2"}, { "P01375", "TNF"}, { "O14931", "NCR3"},
      { "Q9UNN8", "PROCR"}, { "P02730", "SLC4A1"}, { "Q9NSE2", "CISH"}, { "Q96A59", "MARVELD3"}, { "Q9Y231", "FUT9"},
      { "P19320", "VCAM1"}, { "P58753", "TIRAP"}, { "P04921", "GYPC"}, { "P0C091", "FREM3"}, { "P02724", "GYPA"},
      { "P11413", "G6PD"}, { "Q8N126", "CADM3"},  { "Q16570", "ACKR1"}, { "P23634", "ATP2B4"}, { "P17927", "CR1"},
      { "P16442", "ABO"}, {"P69905", "HBA1"}, {"P35613", "BSG"}, {"P08174", "CD55"}, {"Q8NHL6",  "LILRB1"},
      {"Q6GTX8", "LAIR1"} };

  void initializeOntology(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  kol::OntologySetType<std::string> getMalariaGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr, kol::GO::Ontology ontology);
  std::vector<std::string> getAllGOVector(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr, kol::GO::Ontology ontology);
  std::shared_ptr<const kol::LinSimilarity> getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);
  std::shared_ptr<const kol::ResnikSimilarity> getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);

  // Cached GO Term and Gene Data.
  [[nodiscard]] const std::vector<std::string>& malariaGeneVector() const { return malaria_gene_vector_; }
  [[nodiscard]] static const std::map<std::string, std::string>& malariaGeneMap() { return malaria_gene_map_; }

  [[nodiscard]] std::vector<std::string> malariaGOTermsBP() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(malaria_go_terms_BP_.begin(), malaria_go_terms_BP_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsBP() const { return all_go_terms_BP_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache& cacheBP() const { return *cache_BP_ptr_; }

  [[nodiscard]] std::vector<std::string> malariaGOTermsMF() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(malaria_go_terms_MF_.begin(), malaria_go_terms_MF_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsMF() const { return all_go_terms_MF_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache& cacheMF() const { return *cache_MF_ptr_; }

  [[nodiscard]] std::vector<std::string> malariaGOTermsCC() const {

    std::vector<std::string> malaria_GO_vec;
    std::copy(malaria_go_terms_CC_.begin(), malaria_go_terms_CC_.end(), std::back_inserter(malaria_GO_vec));
    return malaria_GO_vec;

  }
  [[nodiscard]] const std::vector<std::string>& allGOTermsCC() const { return all_go_terms_CC_; }
  [[nodiscard]] const kol::AsymmetricSimilarityCache&  cacheCC() const { return *cache_CC_ptr_; }

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class OntologyStats {

public:

  OntologyStats() = default;
  OntologyStats(const OntologyStats &) = default;
  ~OntologyStats() = default;

  OntologyStats &operator=(const OntologyStats &) = default;

  void writeOntology(std::ostream& out_file, char output_delimiter) const;
  void writeOntologyHeader(std::ostream& out_file, char output_delimiter) const;
  void processOntologyStats(const GeneCharacteristic& gene_info,
                            const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                            const OntologyCache& ontology_cache);


private:

  double score_BP_{0.0};
  double score_MF_{0.0};
  double score_CC_{0.0};
  double max_score_{0.0};
  double av_score_{0.0};
  bool targetGene_{false};
};


} // namespace


#endif //KGL_KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
