//
// Created by kellerberrin on 27/3/21.
//

#include "ggtk.hpp"
#include <utility>

#include <boost/filesystem.hpp>


auto test_ggtk(const std::string &goa_file) {

  GoaAnnotationParser annotation_parser;

  return annotation_parser.parseAnnotationFile(goa_file);

}

static const constexpr char* DIRECTORY_ARG = "--dir";
static const constexpr char* GO_OBO = "example_graphs/go-basic.obo";
static const constexpr char* GO_XML = "example_graphs/go_daily-termdb.obo-xml";


int main(int argc, char const **argv) {

  std::string work_directory;
  std::string go_obo = GO_OBO;
  std::string go_xml = GO_XML;

  for (size_t arg = 0; std::cmp_less(arg, argc); ++arg) {

    if (std::string(argv[arg]) == std::string(DIRECTORY_ARG)) {

      ++arg;
      if (std::cmp_less(arg, argc)) {

        work_directory = argv[arg];
        go_obo = work_directory + "/" + go_obo;
        go_xml = work_directory + "/" + go_xml;

      }

    }

  }


  auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
  std::cout << "go_obo file: " << go_obo << ", valid obo file: " << (go_parser_ptr->isFileGood(go_obo) ? "TRUE" : "FALSE") << '\n';
  std::shared_ptr<const GoGraph> obo_go_graph = go_parser_ptr->parseGoFile(go_obo);
  std::cout << "obo file go graph vertices: " << obo_go_graph->getNumVertices() << ", edges: " << obo_go_graph->getNumEdges() << '\n';

  auto xml_parser_ptr = GoParserFactory::createGoParser(GoParserType::XML_GO_STANDARD);
  std::cout << "go_xml file: " << go_xml << ", valid xml file: " << (xml_parser_ptr->isFileGood(go_xml) ? "TRUE" : "FALSE") << '\n';
  std::shared_ptr<const GoGraph> xml_go_graph = xml_parser_ptr->parseGoFile(go_xml);
  std::cout <<  "xml file go graph vertices: " << xml_go_graph->getNumVertices()  << ", edges: " << xml_go_graph->getNumEdges() << '\n';

}


