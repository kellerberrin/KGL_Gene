//
// Created by kellerberrin on 16/5/21.
//

#include "kgl_ontology_database_test.h"
#include "kel_exec_env.h"
#include "kgl_gene_app.h"
#include "kel_utility.h"
#include "kol_InformationContentDAG.h"
#include "kol_InformationContent.h"

#include <fstream>


namespace kgl = kellerberrin::genome;



std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getDAGLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::InformationContentDAG> ic_map_ptr(std::make_shared<kol::InformationContentDAG>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityLin>(ic_map_ptr);

}


std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getDAGResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::InformationContentDAG> ic_map_ptr(std::make_shared<kol::InformationContentDAG>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityResnik>(ic_map_ptr);

}

std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getDAGJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::InformationContentDAG> ic_map_ptr(std::make_shared<kol::InformationContentDAG>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityJIangConrath>(ic_map_ptr);

}


std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {


  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityLin>(ic_map_ptr);

}


std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityResnik>(ic_map_ptr);

}

std::shared_ptr<const kol::SimilarityInterface> kgl::OntologyDatabaseTest::getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::InformationContent> ic_map_ptr(std::make_shared<kol::InformationContent>(ontology_db_ptr->goGraph(), ontology_db_ptr->annotation()));

  return std::make_shared<const kol::SimilarityJIangConrath>(ic_map_ptr);

}


std::shared_ptr<const kol::SetSimilarityInterface> kgl::OntologyDatabaseTest::getBestMatchAverage(const std::shared_ptr<const kol::SimilarityInterface>& similarity_ptr) const {

  return std::make_shared<const kol::SetSimilarityBestMatchAverage>(similarity_ptr);

}

void kgl::OntologyDatabaseTest::performTests() const {

  calcPairs();
//  checkICs();
  calcGenePairs();

}

void kgl::OntologyDatabaseTest::calcPairs() const {

  auto lin_calc_ptr = getDAGLinSimilarity(ontology_db_ptr_);
  double compare = lin_calc_ptr->calculateTermSimilarity("GO:0071312", "GO:0071354");
  ExecEnv::log().info("Check DAG Lin ('GO:0071312', 'GO:0071354') for BP : {}", compare);

  lin_calc_ptr = getLinSimilarity(ontology_db_ptr_);
  compare = lin_calc_ptr->calculateTermSimilarity("GO:0071312", "GO:0071354");
  ExecEnv::log().info("Check Unique Lin ('GO:0071312', 'GO:0071354') for BP : {}", compare);

  auto jiang_calc_ptr = getDAGJiangSimilarity(ontology_db_ptr_);
  compare = jiang_calc_ptr->calculateTermSimilarity("GO:0004022", "GO:0005515");
  ExecEnv::log().info("Check DAG Jiang Conrath ('GO:0004022', 'GO:0005515') for MF : {}", compare);

  jiang_calc_ptr = getJiangSimilarity(ontology_db_ptr_);
  compare = jiang_calc_ptr->calculateTermSimilarity("GO:0004022", "GO:0005515");
  ExecEnv::log().info("Check Unique Jiang Conrath ('GO:0004022', 'GO:0005515') for MF : {}", compare);

  auto resnik_calc_ptr = getDAGResnikSimilarity(ontology_db_ptr_);
  compare = resnik_calc_ptr->calculateTermSimilarity("GO:0004022", "GO:0005515");
  ExecEnv::log().info("Check DAG Resnik ('GO:0004022', 'GO:0005515') for MF : {}", compare);

  resnik_calc_ptr = getResnikSimilarity(ontology_db_ptr_);
  compare = resnik_calc_ptr->calculateTermSimilarity("GO:0004022", "GO:0005515");
  ExecEnv::log().info("Check Unique Resnik ('GO:0004022', 'GO:0005515') for MF : {}", compare);


}


void kgl::OntologyDatabaseTest::checkICs() const {

  std::string file_path = Utility::filePath(IC_FILE_NAME_, GeneExecEnv::getArgs().workDirectory);
  std::ifstream ic_file(file_path);
  std::string line;
  std::map<std::string, double> IC_map;

  while (not std::getline(ic_file, line).eof()) {

    std::vector<std::string> tokens = Utility::char_tokenizer(line, ',');

    if (tokens.size() != 2) {

      ExecEnv::log().error("OntologyDatabaseTest::checkICs; file: {}, line: {}, found token count: {}", file_path, line, tokens.size());
      continue;

    }

    std::string go_term = Utility::trimEndWhiteSpace(tokens[0]);
    std::string ic_text = Utility::trimEndWhiteSpace(tokens[1]);
    try {

      double ic_value = std::stod(ic_text);
      auto result = IC_map.try_emplace(go_term, ic_value);
      if (not result.second) {

        ExecEnv::log().error("OntologyDatabaseTest::checkICs; attempt to add duplicate go term: {}", go_term);

      }

    } catch(...) {

      ExecEnv::log().warn("OntologyDatabaseTest::checkICs; go term: {}, non-numeric ic value text: {}", go_term, ic_text);

    }


  }

  ExecEnv::log().info("Loaded: {} go_term, ic_value pairs from: {}", IC_map.size(), file_path);

  std::string out_file_path =   Utility::filePath(std::string("Out_") + std::string(IC_FILE_NAME_) , GeneExecEnv::getArgs().workDirectory);
  std::ofstream out_file(out_file_path);
  if (not out_file.good()) {

    ExecEnv::log().error("OntologyDatabaseTest::checkICs; problem opening output file: {}", out_file_path);

  }

  std::shared_ptr<const kol::InformationContentDAG> term_map_ptr_ = std::make_shared<const kol::InformationContentDAG>(ontology_db_ptr_->goGraph(),
                                                                                                                       ontology_db_ptr_->annotation());

  std::shared_ptr<const kol::InformationContent> unique_map_ptr_ = std::make_shared<const kol::InformationContent>(ontology_db_ptr_->goGraph(),
                                                                                                                   ontology_db_ptr_->annotation());

  for (auto const& [go_term, ic_value] : IC_map) {

    double alt_ic_value = term_map_ptr_->termInformation(go_term);
    double unique_ic_value = unique_map_ptr_->termInformation(go_term);
    out_file << go_term << ',' << ic_value << ',' << unique_ic_value << ',' << alt_ic_value <<  '\n';

  }

}

void kgl::OntologyDatabaseTest::calcGenePairs() const {

  std::ofstream gene_test_file(GENE_FILE_NAME_);

  if (not gene_test_file.good()) {

    ExecEnv::log().error("OntologyDatabaseTest::calcGenePairs; could not opemn file: {}", GENE_FILE_NAME_);
    return;

  }

  gene_test_file << "gene_row"
                 << FIELD_DELIMITER_
                 << "gene_col"
                 << FIELD_DELIMITER_
                 << "Lin_BP"
                 << FIELD_DELIMITER_
                 << "DAG_Lin_BP"
                 << FIELD_DELIMITER_
                 << "Resnik_MF"
                 << FIELD_DELIMITER_
                 << "DAG_Resnik_MF"
                 << FIELD_DELIMITER_
                 << "Jiang_CC"
                 << FIELD_DELIMITER_
                 << "DAG_Jiang_CC"
                 << '\n';

  std::shared_ptr<const kol::TermAnnotation> term_annotation_ptr(std::make_shared<const kol::TermAnnotation>(ontology_db_ptr_->annotation()->getAllGAFRecords(),
                                                                                                             kol::AnnotationGeneName::SYMBOLIC_GENE_ID));
  auto Lin_BP_ptr = getBestMatchAverage(getLinSimilarity(ontology_db_ptr_));
  auto Resnik_MF_ptr = getBestMatchAverage(getResnikSimilarity(ontology_db_ptr_));
  auto Jiang_CC_ptr = getBestMatchAverage(getJiangSimilarity(ontology_db_ptr_));

  auto DAG_Lin_BP_ptr = getBestMatchAverage(getDAGLinSimilarity(ontology_db_ptr_));
  auto DAG_Resnik_MF_ptr = getBestMatchAverage(getDAGResnikSimilarity(ontology_db_ptr_));
  auto DAG_Jiang_CC_ptr = getBestMatchAverage(getDAGJiangSimilarity(ontology_db_ptr_));

  // Generate a list of genes.
  std::vector<std::string> gene_vector;
  for (auto const & [gaf_if, symbolic_id] : malaria_gene_map_) {

    gene_vector.push_back(symbolic_id);

  }

  for (auto const& gene_row : gene_vector) {

    for (auto const& gene_col : gene_vector) {

      auto gene_row_terms_BP = term_annotation_ptr->getGoTermsForGeneBP(gene_row);
      auto gene_col_terms_BP = term_annotation_ptr->getGoTermsForGeneBP(gene_col);
      double Lin_BP = Lin_BP_ptr->calculateSimilarity(gene_row_terms_BP, gene_col_terms_BP);
      double DAG_Lin_BP = DAG_Lin_BP_ptr->calculateSimilarity(gene_row_terms_BP, gene_col_terms_BP);

      auto gene_row_terms_MF = term_annotation_ptr->getGoTermsForGeneMF(gene_row);
      auto gene_col_terms_MF = term_annotation_ptr->getGoTermsForGeneMF(gene_col);
      double Resnik_MF = Resnik_MF_ptr->calculateSimilarity(gene_row_terms_MF, gene_col_terms_MF);
      double DAG_Resnik_MF = DAG_Resnik_MF_ptr->calculateSimilarity(gene_row_terms_MF, gene_col_terms_MF);

      auto gene_row_terms_CC = term_annotation_ptr->getGoTermsForGeneCC(gene_row);
      auto gene_col_terms_CC = term_annotation_ptr->getGoTermsForGeneCC(gene_col);
      double Jiang_CC = Jiang_CC_ptr->calculateSimilarity(gene_row_terms_CC, gene_col_terms_CC);
      double DAG_Jiang_CC = DAG_Jiang_CC_ptr->calculateSimilarity(gene_row_terms_CC, gene_col_terms_CC);

      gene_test_file << gene_row
                     << FIELD_DELIMITER_
                     << gene_col
                     << FIELD_DELIMITER_
                     << Lin_BP
                     << FIELD_DELIMITER_
                     << DAG_Lin_BP
                     << FIELD_DELIMITER_
                     << Resnik_MF
                     << FIELD_DELIMITER_
                     << DAG_Resnik_MF
                     << FIELD_DELIMITER_
                     << Jiang_CC
                     << FIELD_DELIMITER_
                     << DAG_Jiang_CC
                     << '\n';

    }

  }








}