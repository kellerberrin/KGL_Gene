//
// Created by kellerberrin on 6/4/21.
//


#include "kol_test.h"


struct MICASharedValues {

  const static constexpr double TEST_SHARED_INFORMATION_BAD_IDS{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_DEEP_TERMS_BP{7.919523306};
  const static constexpr double TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM{5.70882396};
  const static constexpr double TEST_SHARED_INFORMATION_CC{5.44229784};
  const static constexpr double TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_CC{14.91515631};
  const static constexpr double TEST_SHARED_INFORMATION_BP{7.077768859};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_BP{15.21457972};
  const static constexpr double TEST_SHARED_INFORMATION_MF{6.50135459};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_MF{12.41485760};
  const static constexpr double TEST_CROSS_ONTOLOGY_SHARED_INFORMATION{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSharedSimilarity<kol::MICASharedInformation, MICASharedValues>;

BOOST_FIXTURE_TEST_SUITE(TestMICASharedSuite, TestSimilarityClass)

#include "kol_test_sharedinfo.h"

BOOST_AUTO_TEST_SUITE_END()
