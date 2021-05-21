//
// Created by kellerberrin on 16/5/21.
//

#ifndef KGL_ONTOLOGY_DATABASE_TEST_H
#define KGL_ONTOLOGY_DATABASE_TEST_H

#include "kol_SimilarityLin.h"
#include "kol_SimilarityResnik.h"
#include "kol_SimilarityJiangConrath.h"
#include "kol_BestMatchAverageSetSimilarity.h"
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

  [[nodiscard]] std::shared_ptr<const kol::LinSimilarity> getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::ResnikSimilarity> getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::JiangConrathSimilarity> getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::LinSimilarity> getLinUnique(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::ResnikSimilarity> getResnikUnique(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;
  [[nodiscard]] std::shared_ptr<const kol::JiangConrathSimilarity> getJiangUnique(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const;

  void calcPairs() const;
  void checkICs() const;

};


} // Namespace.


#endif //KGL_ONTOLOGY_DATABASE_TEST_H
