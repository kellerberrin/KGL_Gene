//
// Created by kellerberrin on 7/4/21.
//
#include "kel_exec_env.h"
#include "kol_test.h"
#include "kol_InformationContent.h"

namespace kel = kellerberrin;
namespace kellerberrin::ontology {

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the various analysis objects as static pointers so they are only created once.
/////////////////////////////////////////////////////////////////////////////////////////
class TestCache_LIN_BP {

public:

  TestCache_LIN_BP() = default;

  ~TestCache_LIN_BP() = default;

  [[nodiscard]] static const std::shared_ptr<const SimilarityLin> &linSimPtr() {

    if (not term_similarity_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr);
    return term_similarity_ptr;

  }

  [[nodiscard]] static const std::shared_ptr<const GoGraphImpl> &goGraphPtr() {

    if (not go_graph_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(go_graph_ptr);
    return go_graph_ptr;

  }

  [[nodiscard]] static const std::shared_ptr<const TermAnnotation> &annotationPtr() {

    if (not annotation_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(annotation_ptr);
    return annotation_ptr;

  }

  [[nodiscard]] static const TermSimilarityCache &cacheSimilar() {

    if (not cache_similarity_ptr) {

      getCacheAnalysis();

    }
    BOOST_REQUIRE(cache_similarity_ptr);
    return *cache_similarity_ptr;

  }

  [[nodiscard]] static const std::vector<std::string>& malariaGenes() { return malaria_genes_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

private:

  inline static std::shared_ptr<const SimilarityLin> term_similarity_ptr;
  inline static std::shared_ptr<const TermSimilarityCache> cache_similarity_ptr;
  inline static std::shared_ptr<const GoGraphImpl> go_graph_ptr;
  inline static std::shared_ptr<const TermAnnotation> annotation_ptr;

  [[nodiscard]] static std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::shared_ptr<const GoGraphImpl> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  static void getCacheAnalysis() {

    go_graph_ptr = getGoGraph();
    annotation_ptr = getAnnotation();
    std::shared_ptr<const InformationContent> ic_map_ptr(std::make_shared<const InformationContent>(go_graph_ptr, annotation_ptr));
    term_similarity_ptr = std::make_shared<const SimilarityLin>(ic_map_ptr);
    cache_similarity_ptr = std::make_shared<const TermSimilarityCache>(annotation_ptr, term_similarity_ptr, GO::Ontology::BIOLOGICAL_PROCESS);

  }

  inline static const std::vector<std::string> malaria_genes_ {
      "P16671", "P06028", "P12318", "P31994", "P05362",
      "O14931", "P68871", "P35228", "P01375", "O14931",
      "Q9UNN8", "P02730", "Q9NSE2", "Q96A59", "Q9Y231",
      "P19320", "P58753", "P04921", "P0C091", "P02724",
      "P11413", "Q8N126", "Q16570", "P23634", "P17927",
      "P16442", "P69905", "P35613", "P08174", "Q8NHL6",
      "Q6GTX8" };


};

} // namespace

namespace kol = kellerberrin::ontology;

BOOST_FIXTURE_TEST_SUITE(TestCache_LIN_Suite, kol::TestCache_LIN_BP)


//////////////////////////////////////////////////////////////////////////
// Test SimilarityMatrix raises error on bad filename
/////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_precomputed_cache_bad_ontology)
{

  auto cache_ptr = std::make_unique<const kol::TermSimilarityCache>(annotationPtr(), linSimPtr(), kol::GO::Ontology::ONTO_ERROR);
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

  auto const ontology_terms = annotationPtr()->getOntologyTerms(kol::GO::Ontology::BIOLOGICAL_PROCESS);
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

BOOST_AUTO_TEST_CASE(test_asymmetric_similarity_Lin_BP)
{

  auto const ontology_terms = annotationPtr()->getOntologyTerms(kol::GO::Ontology::BIOLOGICAL_PROCESS);
  kol::OntologySetType<std::string> go_terms;
  for (auto const& gene : malariaGenes()) {

    go_terms = kol::SetUtilities::setUnion(go_terms, annotationPtr()->getGoTermsForGeneBP(gene));

  }
  std::vector<std::string> malaria_terms = kol::SetUtilities::convertSet(go_terms);
  auto asymmetric_cache_ptr(std::make_shared<const kol::SimilarityCacheAsymmetric>(malaria_terms, ontology_terms, linSimPtr()));
  size_t term_look_up{0};
  for (auto const& row_term : malaria_terms) {

    for (auto const& column_term : ontology_terms) {

      if (cacheSimilar().calculateTermSimilarity(row_term, column_term) != asymmetric_cache_ptr->calculateTermSimilarity(row_term, column_term)) {

        BOOST_FAIL("Term Similarity Cache and Asymmtetric Cache Different row: " + row_term + " column: "
        + column_term + " cache value: " + std::to_string(cacheSimilar().calculateTermSimilarity(row_term, column_term)) +
        " asymmetric value: " + std::to_string(asymmetric_cache_ptr->calculateTermSimilarity(row_term, column_term)) +
        " lin value: " + std::to_string(linSimPtr()->calculateTermSimilarity(row_term, column_term)));

      }

      term_look_up += 2;

    }

  }
  BOOST_TEST_MESSAGE( "test_asymmetric_similarity_Lin_BP ... OK, Looked up terms: " + std::to_string(term_look_up));

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

