//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsMaxSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.53039471027009177};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.53039471027009177};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.51160913084658399};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.65650687706855415};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.65650687706855415};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.57305365849318746};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.48085467210274191};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.48085467210274191};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.40694935907978297};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::SimilarityResnik, kol::SetSimilarityAllPairsMax, AllPairsMaxSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsMaxSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
