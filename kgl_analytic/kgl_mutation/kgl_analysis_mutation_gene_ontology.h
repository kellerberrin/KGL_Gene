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

  OntologyCache( const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
                 const std::shared_ptr<const kol::GoGraph>& go_graph_ptr);
  OntologyCache(const OntologyCache &) = delete;
  ~OntologyCache();

  OntologyCache &operator=(const OntologyCache &) = delete;


  [[nodiscard]] double setSimilarityBP(const std::string& gene) const;
  [[nodiscard]] double setSimilarityMF(const std::string& gene) const;
  [[nodiscard]] double setSimilarityCC(const std::string& gene) const;
  [[nodiscard]] bool isTargetGene(const std::string& gene) const;

private:

  std::shared_ptr<const OntologyGeneCache> gene_cache_ptr_;

  // From the OMIM entry #611162 available at page https://www.omim.org/entry/611162
  inline static const std::map<std::string, std::string> malaria_gene_map_ {
      { "P16671", "CD36"}, { "P06028","GYPB"}, { "P12318", "FCGR2A"}, { "P31994", "FCGR2B"}, { "P05362", "ICAM1"},
      {  "O14931", "NCR3"}, { "P68871", "HBB"}, { "P35228", "NOS2"}, { "P01375", "TNF"}, { "O14931", "NCR3"},
      { "Q9UNN8", "PROCR"}, { "P02730", "SLC4A1"}, { "Q9NSE2", "CISH"}, { "Q96A59", "MARVELD3"}, { "Q9Y231", "FUT9"},
      { "P19320", "VCAM1"}, { "P58753", "TIRAP"}, { "P04921", "GYPC"}, { "P0C091", "FREM3"}, { "P02724", "GYPA"},
      { "P11413", "G6PD"}, { "Q8N126", "CADM3"},  { "Q16570", "ACKR1"}, { "P23634", "ATP2B4"}, { "P17927", "CR1"},
      { "P16442", "ABO"}, {"P69905", "HBA1"}, {"P35613", "BSG"}, {"P08174", "CD55"}, {"Q8NHL6",  "LILRB1"},
      {"Q6GTX8", "LAIR1"} };


  void initializeOntology( const std::shared_ptr<const kol::TermAnnotation>& annotation_ptr,
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

};


} // namespace


#endif //KGL_ANALYSIS_MUTATION_GENE_ONTOLOGY_H
