//
// Created by kellerberrin on 3/4/21.
//

#include "kol_test.h"

struct JiangConrathSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_SIMILARITY_CC{0.57283056264298948};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.6278447852796889};
  const static constexpr double TEST_SIMILARITY_BP{0.62649132831332688};
  const static constexpr double TEST_SIMILARITY_MF{0.48126940869976675};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestTermSimilarity<kol::SimilarityJIangConrath, JiangConrathSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestJiangConrathSimSuite, TestSimilarityClass)

#include "kol_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
