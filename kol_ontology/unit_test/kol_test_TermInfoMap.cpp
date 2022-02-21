//
// Created by kellerberrin on 6/4/21.
//


#include "kol_test.h"
#include "contrib/kol_GoGraphImpl.h"
#include <boost/test/unit_test.hpp>

namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestTermInfoMap {

public:

  TestTermInfoMap() = default;

  ~TestTermInfoMap() = default;

  [[nodiscard]] static const InformationContentDAG &termMap() {

    if (not term_map_ptr_) {

      getTermInfoMap();

    }
    BOOST_REQUIRE(term_map_ptr_);
    return *term_map_ptr_;

  }

  [[nodiscard]] static const GoGraph &goGraph() {

    if (not go_graph_ptr_) {

      getTermInfoMap();

    }
    BOOST_REQUIRE(go_graph_ptr_);
    return *go_graph_ptr_;

  }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

private:

  [[nodiscard]] static std::shared_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }


  static void getTermInfoMap() {

    go_graph_ptr_ = getGoGraph();
    BOOST_REQUIRE(go_graph_ptr_);
    PolicyEvidence default_evidence;
    std::shared_ptr<const TermAnnotation> annotation_ptr = ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());
    term_map_ptr_ = std::make_shared<const InformationContentDAG>(go_graph_ptr_, annotation_ptr);

  }

  inline static std::shared_ptr<const InformationContentDAG> term_map_ptr_;
  inline static std::shared_ptr<const GoGraph> go_graph_ptr_;

};

} // namespace

namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestTermInfoMapSuite, kol::TestTermInfoMap)


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_number_of_values)
{

  const size_t vertex_count = goGraph().getGoGraphImpl().getNumVertices();
  if (termMap().getValues().size() != vertex_count) BOOST_FAIL("Info Map Values not equal to GoGraphImpl Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_values ... OK" );

}

BOOST_AUTO_TEST_CASE(test_number_of_keys)
{

  const size_t vertex_count = goGraph().getGoGraphImpl().getNumVertices();
  if (termMap().getValues().size() != vertex_count) BOOST_FAIL("Info Map Keys not equal to GoGraphImpl Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_keys ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin_bad_id)
{

  if (termMap().termInformation("bad_go_term") != 0.0) BOOST_FAIL("Info Map contains bad term  id" );
  BOOST_TEST_MESSAGE( "ttest_contains_builtin_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin)
{

  if (termMap().termInformation("GO:0000015") == 0.0) BOOST_FAIL("Info Map cannot find term" );
  BOOST_TEST_MESSAGE( "test_contains_builtin ... OK" );

}

BOOST_AUTO_TEST_CASE(test_bracket_access)
{

  const double info{12.907426682529115};
  auto value = termMap().termInformation("GO:0000015");
  BOOST_CHECK_CLOSE( value, info, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_bracket_access ... OK" );

}


BOOST_AUTO_TEST_CASE(test_root_depth_0_BP)
{

  if (termMap().termInformation(kol::GO::getRootTermBP()) != 0.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_MF)
{

  if (termMap().termInformation(kol::GO::getRootTermMF()) != 0.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_CC)
{

  if (termMap().termInformation(kol::GO::getRootTermCC()) != 0.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_CC ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
