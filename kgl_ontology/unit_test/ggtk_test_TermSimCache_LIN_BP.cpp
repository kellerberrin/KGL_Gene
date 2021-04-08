//
// Created by kellerberrin on 7/4/21.
//
#include "ggtk_test.h"

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the various analysis objects as static pointers so they are only created once.
/////////////////////////////////////////////////////////////////////////////////////////
class TestCache_LIN_BP
{

public:

  TestCache_LIN_BP() = default;
  ~TestCache_LIN_BP() = default;

  [[nodiscard]] static const std::shared_ptr<const LinSimilarity>& linSimPtr() {

    if (not term_similarity_ptr_) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return term_similarity_ptr_;

  }

  [[nodiscard]] static const std::shared_ptr<const GoGraph>& goGraphPtr() {

    if (not go_graph_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(go_graph_ptr);
    return go_graph_ptr;

  }

  [[nodiscard]] static const std::shared_ptr<const AnnotationData>& annotationPtr() {

    if (not annotation_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(annotation_ptr);
    return annotation_ptr;

  }

  [[nodiscard]] static const TermSimilarityCache& cacheSimilar() {

    if (not cache_similarity_ptr_) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(cache_similarity_ptr_);
    return *cache_similarity_ptr_;

  }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

private:

  inline static std::shared_ptr<const LinSimilarity> term_similarity_ptr_;
  inline static std::shared_ptr<const TermSimilarityCache> cache_similarity_ptr_;
  inline static std::shared_ptr<const GoGraph> go_graph_ptr;
  inline static std::shared_ptr<const AnnotationData> annotation_ptr;

  [[nodiscard]] static std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  static void getCacheAnalysis() {

    go_graph_ptr = getGoGraph();
    annotation_ptr = getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(go_graph_ptr, annotation_ptr));
    term_similarity_ptr_ = std::make_shared<const LinSimilarity>(go_graph_ptr, info_map_ptr);
    for (size_t count = 0; count < 10; ++count) {

      cache_similarity_ptr_ = std::make_shared<const TermSimilarityCache>(go_graph_ptr, annotation_ptr, term_similarity_ptr_, GO::Ontology::CELLULAR_COMPONENT);
      cache_similarity_ptr_ = std::make_shared<const TermSimilarityCache>(go_graph_ptr, annotation_ptr, term_similarity_ptr_, GO::Ontology::MOLECULAR_FUNCTION);
      cache_similarity_ptr_ = std::make_shared<const TermSimilarityCache>(go_graph_ptr, annotation_ptr, term_similarity_ptr_, GO::Ontology::BIOLOGICAL_PROCESS);

    }

  }

};




BOOST_FIXTURE_TEST_SUITE(TestCache_LIN_Suite, TestCache_LIN_BP)


//////////////////////////////////////////////////////////////////////////
// Test PrecomputedMatrixTermSimilarity raises error on bad filename
/////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_precomputed_cache_bad_ontology)
{

  auto cache_ptr = std::make_unique<const TermSimilarityCache>(goGraphPtr(), annotationPtr(), linSimPtr(), GO::Ontology::ONTO_ERROR);
  BOOST_CHECK_EQUAL(cache_ptr->termCount(), 0.0);
  BOOST_TEST_MESSAGE( "test_precomputed_cache_bad_ontology ... OK" );

}

//////////////////////////////////////////////////////////////////
// Non-existent terms used as input
/////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_Lin_bad_ids)
{

  BOOST_CHECK_EQUAL(linSimPtr()->calculateTermSimilarity("bad_id", "bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(cacheSimilar().calculateTermSimilarity("bad_id","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_Lin_bad_ids ... OK" );

}

BOOST_AUTO_TEST_CASE(test_similarity_Lin_1_bad_1_good_id)
{

  BOOST_CHECK_EQUAL(linSimPtr()->calculateTermSimilarity("GO:0007155", "bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(cacheSimilar().calculateTermSimilarity("GO:0007155","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_Lin_1_bad_1_good_id ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Normalized Similarity [0-1], on terms
///////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_normalized_similarity_Lin_BP_reflexive_sim)
{

  auto const ontology_terms = annotationPtr()->getOntologyTerms(*goGraphPtr(), GO::Ontology::BIOLOGICAL_PROCESS);
  size_t term_look_up{0};
  for (auto const& term_A : ontology_terms) {

    for (auto const& term_B : ontology_terms) {

      if (cacheSimilar().calculateTermSimilarity(term_A, term_B) != cacheSimilar().calculateTermSimilarity(term_B, term_A)) {

        BOOST_FAIL("Term Similarity Cache NOT Reflexive");

      }

      term_look_up += 2;

    }

  }
  BOOST_TEST_MESSAGE( "test_normalized_similarity_Lin_BP_reflexive_sim ... OK, Looked up terms: " + std::to_string(term_look_up));

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_Lin_BP)
{

  double term_value = linSimPtr()->calculateTermSimilarity("GO:0007155", "GO:0044406");
  double cache_value = cacheSimilar().calculateTermSimilarity("GO:0007155", "GO:0044406");
  BOOST_CHECK_CLOSE( term_value, cache_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_Lin_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_normalized_similarity_Lin_BP_1_good_1_root)
{

  double term_value = linSimPtr()->calculateTermSimilarity("GO:0007155", "GO:0008150");
  double cache_value = cacheSimilar().calculateTermSimilarity("GO:0007155", "GO:0008150");
  BOOST_CHECK_CLOSE( term_value, cache_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_Lin_BP_1_good_1_root ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()

