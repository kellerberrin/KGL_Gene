//
// Created by kellerberrin on 6/4/21.
//


#include "kol_test.h"

struct BestPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.87638595311283085};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.64090910308113758};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.93431455927095342};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.60511658596942541};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.95167796503881386};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.83620454724454052};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::SimilarityJIangConrath, kol::SetSimilarityBestMatchAverage, BestPairsAvgSetValues>;

BOOST_FIXTURE_TEST_SUITE(BestPairsAvgSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
