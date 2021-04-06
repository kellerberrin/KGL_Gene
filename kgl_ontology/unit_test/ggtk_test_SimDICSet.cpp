//
// Created by kellerberrin on 6/4/21.
//



#include "ggtk_test.h"

struct SimDICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.639792938};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.448659504};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.94979587};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.51465157};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.933680311};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.528753124};


};

using TestSimilarityClass = TestSetTheorySimilarity<MazanduSimDICSetSimilarity, SimDICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimDICSetSuite, TestSimilarityClass)

#include "ggtk_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
