//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsMaxSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.529139364};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.529139364};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.510371095};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.65664689};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.65664689};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.573227697};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.480874635};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.480874635};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.406859406};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::ResnikSimilarity, kol::AllPairsMaxSetSimilarity, AllPairsMaxSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsMaxSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
