//
// Created by kellerberrin on 5/4/21.
//

#include "kol_test.h"

struct RelevanceSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{0.99669757436887496};
  const static constexpr double TEST_SIMILARITY_CC{0.63212600872438718};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_SIMILARITY_BP{0.71390286143977499};
  const static constexpr double TEST_SIMILARITY_MF{0.67252305245403665};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestTermSimilarity<kol::SimilarityRelevance, RelevanceSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestRelevanceSimSuite, TestSimilarityClass)

#include "kol_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()


