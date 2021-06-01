//
// Created by kellerberrin on 6/4/21.
//

#include "kol_test.h"


struct ExInheritSharedValues {

  const static constexpr double TEST_SHARED_INFORMATION_BAD_IDS{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_DEEP_TERMS_BP{3.5391119850933377};
  const static constexpr double TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM{5.7130980407842458};
  const static constexpr double TEST_SHARED_INFORMATION_CC{2.7234682327945956};
  const static constexpr double TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_CC{14.912074082681379};
  const static constexpr double TEST_SHARED_INFORMATION_BP{3.5372102716694918};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_BP{15.208594815476888};
  const static constexpr double TEST_SHARED_INFORMATION_MF{6.4990002120440069};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_MF{12.409796856084533};
  const static constexpr double TEST_CROSS_ONTOLOGY_SHARED_INFORMATION{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSharedSimilarity<kol::InformationExclusiveInherited, ExInheritSharedValues>;

BOOST_FIXTURE_TEST_SUITE(TestExInheritSharedSuite, TestSimilarityClass)

#include "kol_test_sharedinfo.h"

BOOST_AUTO_TEST_SUITE_END()
