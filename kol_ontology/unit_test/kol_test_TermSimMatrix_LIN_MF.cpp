//
// Created by kellerberrin on 7/4/21.
//

#include "kol_test.h"


namespace kellerberrin::ontology {

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the similarity analysis object in a static pointer so that it is only created once.
/////////////////////////////////////////////////////////////////////////////////////////
class TestMatrix_LIN_MF {

public:

  TestMatrix_LIN_MF() = default;

  virtual ~TestMatrix_LIN_MF() = default;


  [[nodiscard]] const TermSimilarityInterface &termSimilar() {

    if (not term_similarity_ptr_) {

      getMatrixAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] const TermSimilarityInterface &matrixSimilar() {

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

  inline static std::shared_ptr<const TermSimilarityInterface> term_similarity_ptr_;
  inline static std::shared_ptr<const TermSimilarityInterface> matrix_similarity_ptr_;


  [[nodiscard]] std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  void getMatrixAnalysis() {

    std::shared_ptr<const GoGraph> graph_ptr = getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = getAnnotation();
    std::shared_ptr<const TermInformationContentMap> ic_map_ptr(std::make_shared<const TermInformationContentMap>(graph_ptr, annotation_ptr));
    std::shared_ptr<const MICASharedInformation> info_map_ptr(std::make_shared<const MICASharedInformation>( graph_ptr, ic_map_ptr));
    term_similarity_ptr_ = std::make_shared<const LinSimilarity>(info_map_ptr);
    std::shared_ptr<const TermSimilarityWriter> sim_writer_ptr = std::make_shared<const TermSimilarityWriter>(graph_ptr, annotation_ptr);
    if constexpr (WRITER_TEST) {

      sim_writer_ptr->writeSimilarityMatrixMF(term_similarity_ptr_, UnitTestDefinitions::matrixFileNameMF());

    }
    matrix_similarity_ptr_ = std::make_shared<const PrecomputedMatrixTermSimilarity>(UnitTestDefinitions::matrixFileNameMF());

  }

};

} // namespace

namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestMatrix_LIN_MFSuite, kol::TestMatrix_LIN_MF)


//////////////////////////////////////////////////////////////////////////
// Test PrecomputedMatrixTermSimilarity raises error on bad filename
/////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_precomputed_matrix_bad_file)
{

  auto pre_compute_ptr = std::make_unique<const kol::PrecomputedMatrixTermSimilarity>("fake_file.txt");
  BOOST_CHECK_EQUAL(pre_compute_ptr->termCount(), 0.0);
  BOOST_TEST_MESSAGE( "test_precomputed_matrix_bad_file ... OK" );

}

//////////////////////////////////////////////////////////////////
// Non-existent terms used as input
/////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_LinSimilarity_bad_ids)
{

  BOOST_CHECK_EQUAL(termSimilar().calculateTermSimilarity("bad_id","bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(matrixSimilar().calculateTermSimilarity("bad_id","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_LinSimilarity_bad_ids ... OK" );

}

BOOST_AUTO_TEST_CASE(test_similarity_LinSimilarity_1_bad_1_good_id)
{

  BOOST_CHECK_EQUAL(termSimilar().calculateTermSimilarity("GO:0032991","bad_id2"), 0.0);
  BOOST_CHECK_EQUAL(matrixSimilar().calculateTermSimilarity("GO:0032991","bad_id2"), 0.0);
  BOOST_TEST_MESSAGE( "test_similarity_LinSimilarity_1_bad_1_good_id ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Normalized Similarity [0-1], on terms
///////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_normalized_similarity_LinSimilarity_MF_reflexive_sim)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0050662", "GO:0050662");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0050662", "GO:0050662");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_LinSimilarity_MF_reflexive_sim ... OK" );

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_LinSimilarity_MF)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0030170", "GO:0050662");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0030170", "GO:0050662");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_LinSimilarity_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_normalized_similarity_LinSimilarity_MF_1_good_1_root)
{

  double term_value = termSimilar().calculateTermSimilarity("GO:0051192", "GO:0003674");
  double matrix_value = matrixSimilar().calculateTermSimilarity("GO:0051192", "GO:0003674");
  BOOST_CHECK_CLOSE( term_value, matrix_value, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_LinSimilarity_MF_1_good_1_root... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
