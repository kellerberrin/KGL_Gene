//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"


struct AncestorMeanSharedValues {

  const static constexpr double TEST_SHARED_INFORMATION_BAD_IDS{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_DEEP_TERMS_BP{4.2266508352562804};
  const static constexpr double TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM{3.4520569226107551};
  const static constexpr double TEST_SHARED_INFORMATION_CC{2.5183439718404204};
  const static constexpr double TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_CC{14.293721043649006};
  const static constexpr double TEST_SHARED_INFORMATION_BP{3.4820707838137075};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_BP{15.075769645935155};
  const static constexpr double TEST_SHARED_INFORMATION_MF{2.1790065708462216};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_MF{11.975847831148844};
  const static constexpr double TEST_CROSS_ONTOLOGY_SHARED_INFORMATION{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSharedSimilarity<kol::InformationAncestorMean, AncestorMeanSharedValues>;

BOOST_FIXTURE_TEST_SUITE(TestAncestorMeanSharedSuite, TestSimilarityClass)

#include "kol_test_sharedinfo.h"

BOOST_AUTO_TEST_SUITE_END()
