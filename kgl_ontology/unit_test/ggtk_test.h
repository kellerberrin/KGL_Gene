//
// Created by kellerberrin on 2/4/21.
//

#ifndef GGTK_TEST_H
#define GGTK_TEST_H

// Set up the necessary file and directory definitions.
class UnitTestDefinitions {

public:

  // Just static definitions.
  UnitTestDefinitions() = delete;

  [[nodiscard]] static std::string oboFileName() { return std::string(GRAPH_DIRECTORY) + std::string(GO_OBO); }
  [[nodiscard]] static std::string xmlFileName() { return std::string(GRAPH_DIRECTORY) + std::string(GO_XML); }

private:

  static const constexpr char* GRAPH_DIRECTORY = "Additional/ggtk/example_graphs/";
  static const constexpr char* GO_OBO = "go-basic.obo";
  static const constexpr char* GO_XML = "go_daily-termdb.obo-xml";

};

#endif //GGTK_TEST_H
