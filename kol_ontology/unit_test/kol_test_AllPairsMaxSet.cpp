//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"

struct AllPairsMaxSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{0.53333934090867452};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.53333934090867452};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.51264181708134549};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{0.69197341965891856};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.69197341965891856};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.56591824607651742};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{0.4965141303053493};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.4965141303053493};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.38657627968948965};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTermSimilarity<kol::SimilarityResnik, kol::SetSimilarityAllPairsMax, AllPairsMaxSetValues>;

BOOST_FIXTURE_TEST_SUITE(TestAllPairsMaxSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
