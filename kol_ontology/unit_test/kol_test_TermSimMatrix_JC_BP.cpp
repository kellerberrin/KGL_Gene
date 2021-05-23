//
// Created by kellerberrin on 7/4/21.
//
#include "kol_test.h"
#include "kol_InformationContent.h"


namespace kellerberrin::ontology {

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the similarity analysis object in a static pointer so that it is only created once.
/////////////////////////////////////////////////////////////////////////////////////////
class TestMatrix_JC_BP {

public:

  TestMatrix_JC_BP() = default;

  virtual ~TestMatrix_JC_BP() = default;


  [[nodiscard]] const SimilarityInterface &termSimilar() {

    if (not term_similarity_ptr_) {

      getMatrixAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] const SimilarityInterface &matrixSimilar() {

    if (not matrix_similarity_ptr_) {

      getMatrixAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};
  // Enable the test of the matrix writer - this will take 20-25 minutes per test.
  const static constexpr bool WRITER_TEST{false};

private:

  inline static std::shared_ptr<const SimilarityInterface> term_similarity_ptr_;
  inline static std::shared_ptr<const SimilarityInterface> matrix_similarity_ptr_;


  [[nodiscard]] std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = ParserAnnotationFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           PolicyAllowedEvidence());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  void getMatrixAnalysis() {

    std::shared_ptr<const GoGraph> graph_ptr = getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = getAnnotation();
    std::shared_ptr<const InformationContent> ic_map_ptr(std::make_shared<const InformationContent>(graph_ptr, annotation_ptr));
    term_similarity_ptr_ = std::make_shared<const SimilarityJIangConrath>(ic_map_ptr);
    std::shared_ptr<const SimilarityWriter> sim_writer_ptr = std::make_unique<SimilarityWriter>(graph_ptr, annotation_ptr);
    if constexpr (WRITER_TEST) {

      sim_writer_ptr->writeSimilarityMatrixBP(term_similarity_ptr_, UnitTestDefinitions::matrixFileNameBP());

    }
    matrix_similarity_ptr_ = std::make_shared<const SimilarityMatrix>(UnitTestDefinitions::matrixFileNameBP());

  }

};

} // namespace

namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestMatrix_JC_BP_Suite, kol::TestMatrix_JC_BP)


//////////////////////////////////////////////////////////////////////////
// Test SimilarityMatrix raises error on bad filename
/////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_precomputed_matrix_bad_file)
{

  auto pre_compute_ptr = std::make_unique<const kol::SimilarityMatrix>("fake_file.txt");
  BOOST_CHECK_EQUAL(pre_compute_ptr->termCount(), 0.0);
  BOOST_TEST_MESSAGE( "test_precomputed_matrix_bad_file ... OK" );

}

//////////////////////////////////////////////////////////////////
// Non-existent terms used as input
/////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_JiangConrath_bad_ids)
{

  BOOST_CHECK_EQUAL(termSimilar().calculateTermSimilarity("bad_id","bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(matrixSimilar().calculateTermSimilarity("bad_id","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_JiangConrath_bad_ids ... OK" );

}

BOOST_AUTO_TEST_CASE(test_similarity_JiangConrath_1_bad_1_good_id)
{

  BOOST_CHECK_EQUAL(termSimilar().calculateTermSimilarity("GO:0007155","bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(matrixSimilar().calculateTermSimilarity("GO:0007155","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_JiangConrath_1_bad_1_good_id ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Normalized Similarity [0-1], on terms
///////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_normalized_similarity_JiangConrath_BP_reflexive_sim)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0007155", "GO:0007155");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0007155", "GO:0007155");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_JiangConrath_BP_reflexive_sim ... OK" );

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_JiangConrath_BP)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0007155", "GO:0044406");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0007155", "GO:0044406");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_JiangConrath_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_normalized_similarity_JiangConrath_BP_1_good_1_root)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0007155", "GO:0008150");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0007155", "GO:0008150");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_JiangConrath_BP_1_good_1_root ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
