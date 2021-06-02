//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"

struct SimUICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.72323621031173235};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.3106082121140285};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.85361424837233557};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.46338940921070382};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.48746338271902645};


};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityMazanduSimUIC, SimUICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimUICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()

