//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct SimDICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.83939300484046742};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.44816192895599571};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.9210268524013524};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.5023634893294503};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.52406461205320554};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityMazanduSimDIC, SimDICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimDICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
