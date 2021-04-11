//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct SimUISetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.473684210};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.414634146};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.92307692};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.446428571};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.91666666};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.53125};


};


namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestGraphSetSimilarity<kol::GentlemanSimUISetSimilarity, SimUISetValues>;

BOOST_FIXTURE_TEST_SUITE(SimUISetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
