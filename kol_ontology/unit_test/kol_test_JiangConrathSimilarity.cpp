//
// Created by kellerberrin on 3/4/21.
//

#include "kol_test.h"

struct JiangConrathSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_SIMILARITY_CC{0.57983287537685801};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.61688105832177043};
  const static constexpr double TEST_SIMILARITY_BP{0.62827710757723021};
  const static constexpr double TEST_SIMILARITY_MF{0.49232590643658714};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestTermSimilarity<kol::SimilarityJIangConrath, JiangConrathSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestJiangConrathSimSuite, TestSimilarityClass)

#include "kol_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
