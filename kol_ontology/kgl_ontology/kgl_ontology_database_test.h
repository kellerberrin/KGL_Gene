//
// Created by kellerberrin on 16/5/21.
//

#ifndef KGL_ONTOLOGY_DATABASE_TEST_H
#define KGL_ONTOLOGY_DATABASE_TEST_H

#include "kol_SimilarityLin.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityJiangConrath.h"
#include "kol_SetSimilarityBestMatchAverage.h"
#include "kgl_ontology_database.h"

namespace kol = kellerberrin::ontology;
namespace kellerberrin::genome {   //  organization::project level namespace



class OntologyDatabaseTest {

public:

  explicit OntologyDatabaseTest(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) : ontology_db_ptr_(ontology_db_ptr) {}
  ~OntologyDatabaseTest() = default;

  void performTests() const;

private:

  std::shared_ptr<const kol::OntologyDatabase> ontology_db_ptr_;
  const static constexpr char* IC_FILE_NAME_ = "Test_IC.csv";
  const static constexpr char* GENE_FILE_NAME_ = "Test_Gene.csv";
  const static constexpr char FIELD_DELIMITER_ = ',';

  // From the OMIM entry #611162 available at page https://www.omim.org/entry/611162
  inline static const std::map<std::string, std::string> malaria_gene_map_ {
      { "P16671", "CD36"}, { "P06028","GYPB"}, { "P12318", "FCGR2A"}, { "P31994", "FCGR2B"}, { "P05362", "ICAM1"},
      {  "O14931", "NCR3"}, { "P68871", "HBB"}, { "P35228", "NOS2"}, { "P01375", "TNF"}, { "O14931", "NCR3"},
      { "Q9UNN8", "PROCR"}, { "P02730", "SLC4A1"}, { "Q9NSE2", "CISH"}, { "Q96A59", "MARVELD3"}, { "Q9Y231", "FUT9"},
      { "P19320", "VCAM1"}, { "P58753", "TIRAP"}, { "P04921", "GYPC"}, { "P0C091", "FREM3"}, { "P02724", "GYPA"},
      { "P11413", "G6PD"}, { "Q8N126", "CADM3"},  { "Q16570", "ACKR1"}, { "P23634", "ATP2B4"}, { "P17927", "CR1"},
      { "P16442", "ABO"}, {"P69905", "HBA1"}, {"P35613", "BSG"}, {"P08174", "CD55"}, {"Q8NHL6",  "LILRB1"},
      {"Q6GTX8", "LAIR1"} };

  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getDAGLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getDAGResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getDAGJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SimilarityInterface> getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SetSimilarityInterface> getBestMatchAverage(const std::shared_ptr<const kol::SimilarityInterface>& similarity_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::SetSimilarityInterface> getBMA_Alt(const std::shared_ptr<const kol::SimilarityInterface>& similarity_ptr) const;

  void calcPairs() const;
  void checkICs() const;
  void calcGenePairs() const;

};


} // Namespace.


#endif //KGL_ONTOLOGY_DATABASE_TEST_H
