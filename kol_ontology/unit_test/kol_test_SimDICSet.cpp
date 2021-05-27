//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct SimDICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.64027054299215747};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.44890534304156299};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.94977800963399928};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.51447502071732909};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.93367346653933092};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.52872540686382308};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityMazanduSimDIC, SimDICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimDICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
