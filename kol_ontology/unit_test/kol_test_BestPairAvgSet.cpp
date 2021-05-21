//
// Created by kellerberrin on 6/4/21.
//


#include "kol_test.h"

struct BestPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.90088932933096766};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.64187989536615808};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.93194665064307958};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.59221890479635819};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.97004529765823833};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.83845710341523683};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::JiangConrathSimilarity, kol::BestMatchAverageSetSimilarity, BestPairsAvgSetValues>;

BOOST_FIXTURE_TEST_SUITE(BestPairsAvgSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
