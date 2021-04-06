//
// Created by kellerberrin on 6/4/21.
//




#include "ggtk_test.h"

struct SimGICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.470364371};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.289207627};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.90439167};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.34648541};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.87561012};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.35939116};


};

using TestSimilarityClass = TestSetTheorySimilarity<PesquitaSimGICSetSimilarity, SimGICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimGICSetSuite, TestSimilarityClass)

#include "ggtk_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
