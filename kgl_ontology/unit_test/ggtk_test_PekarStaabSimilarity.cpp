//
// Created by kellerberrin on 5/4/21.
//


#include "ggtk_test.h"



struct PekarStaabSimValues {

  const static constexpr double TEST_SIMILARITY_BAD_IDS{0.0};
  const static constexpr double TEST_SIMILARITY_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_SIMILARITY_CC{0.25};
  const static constexpr double TEST_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_REFLEXIVE_SIM{1.0};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC{0.25};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_SIMILARITY_BP{0.3333333};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_BP{0.3333333};
  const static constexpr double TEST_SIMILARITY_MF{0.5};
  const static constexpr double TEST_NORMALIZED_SIMILARITY_MF{0.5};
  const static constexpr double TEST_CROSS_ONTOLOGY_SIMILARITY{0.0};
  const static constexpr double TEST_CROSS_ONTOLOGY_NORMALIZED_SIMILARITY{0.0};

};

using TestSimilarityClass = TestPekarStaabSimilarity<PekarStaabSimValues>;

BOOST_FIXTURE_TEST_SUITE(TestPekarStaabSimSuite, TestSimilarityClass)

#include "ggtk_test_similarity.h"

BOOST_AUTO_TEST_SUITE_END()
