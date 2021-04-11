//
// Created by kellerberrin on 5/4/21.
//


#include "kol_test.h"


struct ResnikSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{5.70882396};
  const static constexpr double TEST_SIMILARITY_CC{5.44229784};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_REFLEXIVE_SIM{0.382753210};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC{0.364883728};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_SIMILARITY_BP{7.07776885};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_BP{0.4651964752};
  const static constexpr double TEST_SIMILARITY_MF{6.50135459};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_MF{0.523675325};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};
  const static constexpr double TEST_CROSS_ONTOLOGY_NORMALIZED_SIMILARITY{0.0};

};


namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestTermSimilarity<kol::ResnikSimilarity, ResnikSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestResnikSimilaritySuite, TestSimilarityClass)

#include "kol_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
