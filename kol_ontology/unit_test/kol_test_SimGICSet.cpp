//
// Created by kellerberrin on 6/4/21.
//




#include "kol_test.h"

struct SimGICSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.72323621031173235};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.28879426102395678};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.85361424837233557};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.33543752823208284};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.35507286859097748};


};


namespace kol = kellerberrin::ontology;


using TestSimilarityClass = kol::TestSetTheorySimilarity<kol::SetSimilarityPesquitaSimGIC, SimGICSetValues>;


BOOST_FIXTURE_TEST_SUITE(SimGICSetSuite, TestSimilarityClass)

#include "kol_test_setsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
