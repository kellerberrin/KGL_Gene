//
// Created by kellerberrin on 3/4/21.
//

#include "ggtk_test.h"

struct JiangConrathSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_SIMILARITY_CC{0.789688805};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.808623394};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC{0.789688805};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_1_GOOD_1_ROOT{0.808623394};
  const static constexpr double TEST_SIMILARITY_BP{0.814125375};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_BP{0.81412537};
  const static constexpr double TEST_SIMILARITY_MF{0.74620920};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_MF{0.74620920};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};
  const static constexpr double TEST_CROSS_ONTOLOGY_NORMALIZED_SIMILARITY{0.0};

};

using TestSimilarityClass = TestTermSimilarity<JiangConrathSimilarity, JiangConrathSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestJiangConrathSimSuite, TestSimilarityClass)

#include "ggtk_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
