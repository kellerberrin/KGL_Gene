//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.545562458};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.536042188};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.259268405};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.20152102};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.195834099};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.12567561};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.460019367};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.469187743};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.362365916};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::LinSimilarity, AllPairsAverageSetSimilarity, AllPairsAvgSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsAvgSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
