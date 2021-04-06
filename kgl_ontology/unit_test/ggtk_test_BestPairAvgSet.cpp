//
// Created by kellerberrin on 6/4/21.
//


#include "ggtk_test.h"

struct BestPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.950444664};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.820939947};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.9659733253};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.794793986};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.985022648};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.919228551};


};

using TestSimilarityClass = TestSetTermSimilarity<JiangConrathSimilarity, BestMatchAverageSetSimilarity, BestPairsAvgSetValues>;


BOOST_FIXTURE_TEST_SUITE(BestPairsAvgSetSuite, TestSimilarityClass)

#include "ggtk_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
