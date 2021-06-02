//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.54851899351361788};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.56595997595536263};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.26062078744966805};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.22569510936855908};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.23750079790893752};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.17049893747062};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.45453078469398234};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.52270771292530982};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.36121469690765245};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::SimilarityLin, kol::SetSimilarityAllPairsAverage, AllPairsAvgSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsAvgSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
