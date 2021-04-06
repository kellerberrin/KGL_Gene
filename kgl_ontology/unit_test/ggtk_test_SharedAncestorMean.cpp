//
// Created by kellerberrin on 6/4/21.
//



#include "ggtk_test.h"


struct AncestorMeanSharedValues {

  const static constexpr double TEST_SHARED_INFORMATION_BAD_IDS{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_DEEP_TERMS_BP{4.147530971};
  const static constexpr double TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM{3.717040601};
  const static constexpr double TEST_SHARED_INFORMATION_CC{2.72114892};
  const static constexpr double TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_CC{14.9151563};
  const static constexpr double TEST_SHARED_INFORMATION_BP{3.538884429};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_BP{15.21457972};
  const static constexpr double TEST_SHARED_INFORMATION_MF{2.28338684};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_MF{12.41485760};
  const static constexpr double TEST_CROSS_ONTOLOGY_SHARED_INFORMATION{0.0};

};

using TestSimilarityClass = TestSharedSimilarity<AncestorMeanSharedInformation, AncestorMeanSharedValues>;

BOOST_FIXTURE_TEST_SUITE(TestAncestorMeanSharedSuite, TestSimilarityClass)

#include "ggtk_test_sharedinfo.h"

BOOST_AUTO_TEST_SUITE_END()
