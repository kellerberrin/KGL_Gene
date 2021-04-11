//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"
#include <boost/test/unit_test.hpp>


namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestTermProbMap {

public:

  TestTermProbMap() = default;

  ~TestTermProbMap() = default;

  [[nodiscard]] static const TermProbabilityMap &probMap() {

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

  [[nodiscard]] static std::unique_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }


  static void getTermInfoMap() {

    go_graph_ptr_ = getGoGraph();
    BOOST_REQUIRE(go_graph_ptr_);
    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    std::shared_ptr<const AnnotationData> annotation_ptr = anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());
    term_map_ptr_ = std::make_unique<const TermProbabilityMap>(go_graph_ptr_, annotation_ptr);

  }

  inline static std::unique_ptr<const TermProbabilityMap> term_map_ptr_;
  inline static std::shared_ptr<const GoGraph> go_graph_ptr_;

};

} // namespace


namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestTermProbMapSuite, kol::TestTermProbMap)


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_get_values_type)
{

  auto values = probMap().getValues();
  if constexpr(std::is_same<decltype(values), std::vector<double>>::value) {

    BOOST_TEST_MESSAGE( "test_get_values_type ... OK" );

  } else {

    BOOST_FAIL("Expected returned type 'std::vector<double>>'" );

  }

}


BOOST_AUTO_TEST_CASE(test_get_keys_type)
{

  auto keys = probMap().getKeys();
  if constexpr(std::is_same<decltype(keys), std::vector<std::string>>::value) {

    BOOST_TEST_MESSAGE( "test_get_keys_type ... OK" );

  } else {

    BOOST_FAIL("Expected returned type 'std::vector<std::string>'" );

  }

}

BOOST_AUTO_TEST_CASE(test_number_of_values)
{

  const size_t vertex_count = goGraph().getNumVertices();
  if (probMap().getValues().size() != vertex_count) BOOST_FAIL("Info Map Values not equal to GoGraph Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_values ... OK" );

}

BOOST_AUTO_TEST_CASE(test_number_of_keys)
{

  const size_t vertex_count = goGraph().getNumVertices();
  if (probMap().getKeys().size() != vertex_count) BOOST_FAIL("Info Map Keys not equal to GoGraph Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_keys ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin_bad_id)
{

  if (probMap().hasTerm("bad_go_term")) BOOST_FAIL("Info Map contains bad term  id" );
  BOOST_TEST_MESSAGE( "ttest_contains_builtin_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin)
{

  if (not probMap().hasTerm("GO:0000015")) BOOST_FAIL("Info Map cannot find term" );
  BOOST_TEST_MESSAGE( "test_contains_builtin ... OK" );

}

BOOST_AUTO_TEST_CASE(test_bracket_access)
{
  // This probability was given in the python test rig.
  // const double prob_original{1.33117e-06};
  // However, this is the probability as returned by the unmodified ggtk code.
  const double prob{1.3319560907355129e-06};
  auto value = probMap()["GO:0000015"];
  BOOST_CHECK_CLOSE( value, prob, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_bracket_access ... OK" );

}


BOOST_AUTO_TEST_CASE(test_root_depth_0_BP)
{

  if (probMap()[GO::getRootTermBP()] != 1.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_MF)
{

  if (probMap()[GO::getRootTermMF()] != 1.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_CC)
{

  if (probMap()[GO::getRootTermCC()] != 1.0) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_CC ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
