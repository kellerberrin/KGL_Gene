//
// Created by kellerberrin on 3/4/21.
//

#include <ggtk.h>
#include "ggtk_test.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestJiangConrathSim
{

public:

  TestJiangConrathSim() = default;
  ~TestJiangConrathSim() = default;


  [[nodiscard]] static const JiangConrathSimilarity& jiangConrath() {

    if (not jiang_conrath_ptr_) {

      jiang_conrath_ptr_ = getJiangConrath();

    }
    BOOST_REQUIRE(jiang_conrath_ptr_);
    return *jiang_conrath_ptr_;

  }


private:

  inline static std::unique_ptr<const JiangConrathSimilarity> jiang_conrath_ptr_;

  [[nodiscard]] static std::unique_ptr<AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::unique_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  [[nodiscard]] static std::unique_ptr<const JiangConrathSimilarity> getJiangConrath() {

    std::shared_ptr<const GoGraph> graph_ptr = getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(graph_ptr, annotation_ptr));
    return std::make_unique<JiangConrathSimilarity>(JiangConrathSimilarity(graph_ptr, info_map_ptr));

  }

};


BOOST_FIXTURE_TEST_SUITE(TestJiangConrathSimSuite, TestJiangConrathSim)

///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_bad_ids)
{

  if (jiangConrath().calculateTermSimilarity("bad_id","bad_id2") !=0 ) BOOST_FAIL("Non-zero value on bad id");
  BOOST_TEST_MESSAGE( "test_similarity_bad_ids ... OK" );

}


BOOST_AUTO_TEST_CASE(test_similarity_1_bad_1_good_id)
{

  if (jiangConrath().calculateTermSimilarity("GO:0032991","bad_id2") !=0 ) BOOST_FAIL("Non-zero value on bad id");
  BOOST_TEST_MESSAGE( "test_similarity_1_bad_1_good_id ... OK" );

}


////////////////////////////////////////////////////////
// Similarity on CC terms
///////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_CC_reflexive_sim)
{

  auto value = jiangConrath().calculateTermSimilarity("GO:0043234", "GO:0043234");
  BOOST_CHECK_CLOSE( value, 1.0, 0.0001);
  BOOST_TEST_MESSAGE( "test_similarity_CC_reflexive_sim ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
