//
// Created by kellerberrin on 6/4/21.
//




#include "kol_test.h"

struct SimGICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.47088083566352013};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.28941195885609439};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.90435928627147066};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.34632539196059714};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.87559808111421167};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.35936555237917162};


};


namespace kol = kellerberrin::ontology;


using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityPesquitaSimGIC, SimGICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimGICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
