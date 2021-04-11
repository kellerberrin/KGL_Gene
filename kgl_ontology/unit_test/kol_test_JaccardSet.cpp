//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct JaccardSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.6666666666};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.1428571428};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.8333333333};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.153846153};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.75};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.4};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestJaccardSetSimilarity<JaccardSetValues>;

BOOST_FIXTURE_TEST_SUITE(JaccardSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
