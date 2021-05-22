//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"


struct LinSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_SIMILARITY_CC{0.634365349};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC{0.634365349};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_SIMILARITY_BP{0.714509483};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_BP{0.714509483};
  const static constexpr double TEST_SIMILARITY_MF{0.67356674};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_MF{0.67356674};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};
  const static constexpr double TEST_CROSS_ONTOLOGY_NORMALIZED_SIMILARITY{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestTermSimilarity<kol::SimilarityLin, LinSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestLinSimSuite, TestSimilarityClass)

#include "kol_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
