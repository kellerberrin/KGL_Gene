//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsAvgSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.54562114677762175};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.53605748798021413};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.25922484248779792};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.20136259156248076};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.19569972662593485};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.12546134352547408};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.46006026703081754};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.46922105985171819};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.36229399778512517};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::SimilarityLin, kol::SetSimilarityAllPairsAverage, AllPairsAvgSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsAvgSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
