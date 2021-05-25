//
// Created by kellerberrin on 6/4/21.
//

#include "kol_test.h"
#include <boost/test/unit_test.hpp>


namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestDepthMap {

public:

  TestDepthMap() = default;

  ~TestDepthMap() = default;

  [[nodiscard]] static const InformationDepthMap &depthMap() {

    if (not depth_map_ptr_) {

      getDepthMap();

    }
    BOOST_REQUIRE(depth_map_ptr_);
    return *depth_map_ptr_;

  }

  [[nodiscard]] static const GoGraph &goGraph() {

    if (not go_graph_ptr_) {

      getDepthMap();

    }
    BOOST_REQUIRE(go_graph_ptr_);
    return *go_graph_ptr_;

  }

private:

  [[nodiscard]] static std::shared_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }


  static void getDepthMap() {

    go_graph_ptr_ = getGoGraph();
    BOOST_REQUIRE(go_graph_ptr_);
    depth_map_ptr_ = std::make_unique<const InformationDepthMap>(*go_graph_ptr_);

  }

  inline static std::unique_ptr<const InformationDepthMap> depth_map_ptr_;
  inline static std::shared_ptr<const GoGraph> go_graph_ptr_;

};

} // namespace

namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestDepthMapSuite, kol::TestDepthMap)


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_get_values_type)
{

  auto values = depthMap().getValues();
  if constexpr(std::is_same<decltype(values), std::vector<size_t>>::value) {

    BOOST_TEST_MESSAGE( "test_get_values_type ... OK" );

  } else {

    BOOST_FAIL("Expected returned type 'std::vector<TermDepthType>>'" );

  }

}


BOOST_AUTO_TEST_CASE(test_get_keys_type)
{

  auto keys = depthMap().getKeys();
  if constexpr(std::is_same<decltype(keys), std::vector<std::string>>::value) {

    BOOST_TEST_MESSAGE( "test_get_keys_type ... OK" );

  } else {

    BOOST_FAIL("Expected returned type 'std::vector<std::string>'" );

  }

}

BOOST_AUTO_TEST_CASE(test_number_of_values)
{

  const size_t vertex_count = goGraph().getNumVertices();
  if (depthMap().getValues().size() != vertex_count) BOOST_FAIL("Depth Map Values not equal to GoGraph Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_values ... OK" );

}

BOOST_AUTO_TEST_CASE(test_number_of_keys)
{

  const size_t vertex_count = goGraph().getNumVertices();
  if (depthMap().getKeys().size() != vertex_count) BOOST_FAIL("Depth Map Keys not equal to GoGraph Vertices" );
  BOOST_TEST_MESSAGE( "test_number_of_keys ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin_bad_id)
{

  if (depthMap().hasTerm("bad_go_term")) BOOST_FAIL("Depth Map contains bad term  id" );
  BOOST_TEST_MESSAGE( "ttest_contains_builtin_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_contains_builtin)
{

  if (not depthMap().hasTerm("GO:0000015")) BOOST_FAIL("Depth Map cannot find term" );
  BOOST_TEST_MESSAGE( "test_contains_builtin ... OK" );

}

BOOST_AUTO_TEST_CASE(test_bracket_access)
{

  const size_t depth{4};
  if (depthMap().getValue("GO:0000015") != depth) BOOST_FAIL("Depth Map wrong depth returned" );
  BOOST_TEST_MESSAGE( "test_bracket_access ... OK" );

}


BOOST_AUTO_TEST_CASE(test_bracket_access_2)
{

  const size_t depth{3};
  if (depthMap().getValue("GO:0098900") != depth) BOOST_FAIL("Depth Map wrong depth returned" );
  BOOST_TEST_MESSAGE( "test_bracket_access_2 ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_BP)
{

  const size_t depth{0};
  if (depthMap().getValue(kol::GO::getRootTermBP()) != depth) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_MF)
{

  const size_t depth{0};
  if (depthMap().getValue(kol::GO::getRootTermMF()) != depth) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_depth_0_CC)
{

  const size_t depth{0};
  if (depthMap().getValue(kol::GO::getRootTermCC()) != depth) BOOST_FAIL("Root node depth is not zero" );
  BOOST_TEST_MESSAGE( "test_root_depth_0_CC ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
