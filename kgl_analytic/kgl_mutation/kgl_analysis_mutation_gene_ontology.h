//
// Created by kellerberrin on 14/4/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
#define KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H

#include <fstream>
#include "kgl_ontology_database.h"
#include "kgl_analysis_mutation_gene_stats.h"



namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Forward decl of ontology functionality for PIMPL pattern.
class OntologyGeneCache;

class OntologyCache {

public:

  OntologyCache( const std::vector<std::string>& gene_vector,
                 const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                 const std::shared_ptr<const kol::GoGraph>& go_graph_ptr);
  OntologyCache(const OntologyCache &) = delete;
  ~OntologyCache();

  OntologyCache &operator=(const OntologyCache &) = delete;


  [[nodiscard]] double setSimilarityBP(const std::string& gene) const;
  [[nodiscard]] double setSimilarityMF(const std::string& gene) const;
  [[nodiscard]] double setSimilarityCC(const std::string& gene) const;

  [[nodiscard]] bool isTargetGene(const std::string& gene) const;

  [[nodiscard]] std::pair<std::string, double> maxMFSim(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxBPSim(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxCCSim(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxFunMFBPSim(const std::string& gene) const;

  [[nodiscard]] std::pair<std::string, double> maxMFInfo(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxBPInfo(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxCCInfo(const std::string& gene) const;
  [[nodiscard]] std::pair<std::string, double> maxFunMFBPInfo(const std::string& gene) const;

private:

  std::shared_ptr<const OntologyGeneCache> gene_cache_ptr_;


  void initializeOntology( const std::vector<std::string>& gene_vector,
                           const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                           const std::shared_ptr<const kol::GoGraph>& go_graph_ptr);

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
  void processOntologyStats(const std::string& gene_info, const OntologyCache& ontology_cache);

private:

  double score_BP_{0.0};
  double score_MF_{0.0};
  double score_CC_{0.0};
  double max_score_{0.0};
  double av_score_{0.0};

  bool targetGene_{false};
  // Similarity based statistics
  double max_MF_sim_{0.0};
  std::string max_MF_gene_sim_;
  double max_BP_sim_{0.0};
  std::string max_BP_gene_sim_;
  double max_CC_sim_{0.0};
  std::string max_CC_gene_sim_;
  double max_FunMFBP_sim_{0.0};
  std::string max_FunMFBP_gene_sim_;
  // Information base statistics
  double max_MF_info_{0.0};
  std::string max_MF_gene_info_;
  double max_BP_info_{0.0};
  std::string max_BP_gene_info_;
  double max_CC_info_{0.0};
  std::string max_CC_gene_info_;
  double max_FunMFBP_info_{0.0};
  std::string max_FunMFBP_gene_info_;

};


} // namespace


#endif //KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
