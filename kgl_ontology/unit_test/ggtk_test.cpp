//
// Created by kellerberrin on 1/4/21.
//

#include "ggtk.h"
#include <utility>

//#define BOOST_TEST_MAIN  // in only one cpp file
#define BOOST_TEST_MODULE example
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
#include <boost/bind.hpp>
using namespace boost::unit_test;


static const constexpr char* DIRECTORY = "Additional/ggtk/";
static const constexpr char* GO_OBO = "example_graphs/go-basic.obo";
static const constexpr char* GO_XML = "example_graphs/go_daily-termdb.obo-xml";

class TestGoParsers
{

public:

  TestGoParsers()
  {

    //    BOOST_TEST_MESSAGE("setup TestGoParsers");

    go_parser_ptr_ = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr_);
    obo_go_graph_ = go_parser_ptr_->parseGoFile(go_obo);
    BOOST_REQUIRE(obo_go_graph_);

    xml_parser_ptr_ = GoParserFactory::createGoParser(GoParserType::XML_GO_STANDARD);
    BOOST_REQUIRE(xml_parser_ptr_);
    xml_go_graph_ = xml_parser_ptr_->parseGoFile(go_xml);
    BOOST_REQUIRE(xml_go_graph_);

  }
  ~TestGoParsers() = default;

  [[nodiscard]] bool checkObOFile() const {  return go_parser_ptr_->isFileGood(go_obo); }
  [[nodiscard]] size_t checkObOVertices() const {  return obo_go_graph_->getNumVertices(); }
  [[nodiscard]] size_t checkObOEdges() const {  return obo_go_graph_->getNumEdges(); }

  [[nodiscard]] bool checkXmlFile() const {  return xml_parser_ptr_->isFileGood(go_xml); }
  [[nodiscard]] size_t checkXmlVertices() const {  return xml_go_graph_->getNumVertices(); }
  [[nodiscard]] size_t checkXmlEdges() const {  return xml_go_graph_->getNumEdges(); }

private:

  const std::string go_obo{std::string(DIRECTORY) + std::string(GO_OBO)};
  const std::string go_xml{std::string(DIRECTORY) + std::string(GO_XML)};

  std::unique_ptr<const GoParserInterface> go_parser_ptr_;
  std::unique_ptr<const GoParserInterface> xml_parser_ptr_;

  std::shared_ptr<const GoGraph> obo_go_graph_;
  std::shared_ptr<const GoGraph> xml_go_graph_;

};

BOOST_FIXTURE_TEST_SUITE(TestGoParserSuite, TestGoParsers)

BOOST_AUTO_TEST_CASE( test_case1 )
{

  BOOST_WARN( sizeof(int) < 4 );
}

BOOST_AUTO_TEST_CASE( test_case2 )
{
  BOOST_REQUIRE_EQUAL( 1, 2 );
  BOOST_FAIL( "Should never reach this line" );
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE( test_suite2 )

BOOST_AUTO_TEST_CASE( test_case3 )
{
  BOOST_CHECK( true );
}

BOOST_AUTO_TEST_CASE( test_case4 )
{
  BOOST_CHECK( false );
}

BOOST_AUTO_TEST_SUITE_END()



/*
//____________________________________________________________________________//

class test_class {
public:
  void test_method1()
  {
    BOOST_CHECK( false );
  }
  void test_method2()
  {
    BOOST_CHECK( false );
  }
};

//____________________________________________________________________________//

test_suite*
init_unit_test_suite( int argc, char* argv[] )
{


  std::shared_ptr<test_class> tester(std::make_shared<test_class>());

  framework::master_test_suite().
      add( BOOST_TEST_CASE( [tester] { tester->test_method1(); }));

  framework::master_test_suite().
      add( BOOST_TEST_CASE( [tester] { tester->test_method2(); }));

  return 0;

}

*/

/*

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

*/

