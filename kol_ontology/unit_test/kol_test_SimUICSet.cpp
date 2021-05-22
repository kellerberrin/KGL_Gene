//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct SimUICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.47036437};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.31107195};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.904391678};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.473657694};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.875610120};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.492563893};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityMazanduSimUIC, SimUICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimUICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()

