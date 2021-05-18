//
// Created by kellerberrin on 16/5/21.
//

#include "kgl_ontology_database_test.h"
#include "kel_exec_env.h"
#include "kgl_gene_app.h"
#include "kel_utility.h"
#include "kol_TermInformationContentMap.h"


namespace kgl = kellerberrin::genome;


std::shared_ptr<const kol::LinSimilarity> kgl::OntologyDatabaseTest::getLinSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::TermInformationContentMap> info_map_ptr(std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr->goGraph(),
                                                                                                                            ontology_db_ptr->annotation()));
  return std::make_shared<const kol::LinSimilarity>(ontology_db_ptr->goGraph(), info_map_ptr);

}


std::shared_ptr<const kol::ResnikSimilarity> kgl::OntologyDatabaseTest::getResnikSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::TermInformationContentMap> info_map_ptr(std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr->goGraph(),
                                                                                                                            ontology_db_ptr->annotation()));
  return std::make_shared<const kol::ResnikSimilarity>(ontology_db_ptr->goGraph(), info_map_ptr);

}

std::shared_ptr<const kol::JiangConrathSimilarity> kgl::OntologyDatabaseTest::getJiangSimilarity(const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr) const {

  std::shared_ptr<const kol::TermInformationContentMap> info_map_ptr(std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr->goGraph(),
                                                                                                                            ontology_db_ptr->annotation()));
  return std::make_shared<const kol::JiangConrathSimilarity>(ontology_db_ptr->goGraph(), info_map_ptr);

}


void kgl::OntologyDatabaseTest::performTests() const {

  calcPairs();
  checkICs();

}

void kgl::OntologyDatabaseTest::calcPairs() const {

  auto lin_calc_ptr = getLinSimilarity(ontology_db_ptr_);
  double compare = lin_calc_ptr->calculateTermSimilarity("GO:0071312", "GO:0071354");
  ExecEnv::log().info("Check Lin ('GO:0071312', 'GO:0071354') for BP : {}", compare);
  auto jiang_calc_ptr = getJiangSimilarity(ontology_db_ptr_);
  compare = jiang_calc_ptr->calculateTermSimilarity("GO:0004022", "GO:0005515");
  ExecEnv::log().info("Check Jiang Conrath ('GO:0004022', 'GO:0005515') for MF : {}", compare);


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

    ExecEnv::log().error("OntologyDatabaseTest::checkICs; pronlem opening output file: {}", out_file_path);

  }

  std::shared_ptr<const kol::TermInformationContentMap> term_map_ptr_ = std::make_shared<const kol::TermInformationContentMap>(ontology_db_ptr_->goGraph(),
                                                                                                                               ontology_db_ptr_->annotation());

  std::shared_ptr<const kol::TermInformationContentUnique> unique_map_ptr_ = std::make_shared<const kol::TermInformationContentUnique>(ontology_db_ptr_->goGraph(),
                                                                                                                                       ontology_db_ptr_->annotation());

  for (auto const& [go_term, ic_value] : IC_map) {

    double alt_ic_value = term_map_ptr_->getValue(go_term);
    double unique_ic_value = unique_map_ptr_->getValue(go_term);
    out_file << go_term << ',' << ic_value << ',' << unique_ic_value << ',' << alt_ic_value <<  '\n';

  }

}
