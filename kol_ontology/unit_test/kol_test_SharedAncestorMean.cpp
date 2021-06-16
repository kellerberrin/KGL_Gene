//
// Created by kellerberrin on 6/4/21.
//



#include "kol_test.h"


struct AncestorMeanSharedValues {

  const static constexpr double TEST_SHARED_INFORMATION_BAD_IDS{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID{0.0};
  const static constexpr double TEST_SHARED_INFORMATION_DEEP_TERMS_BP{2.3234463927073175};
  const static constexpr double TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM{1.5993399369487016};
  const static constexpr double TEST_SHARED_INFORMATION_CC{1.1288247867592096};
  const static constexpr double TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT{0.0};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_CC{11.135449967400994};
  const static constexpr double TEST_SHARED_INFORMATION_BP{2.1463018481343532};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_BP{11.714379489913354};
  const static constexpr double TEST_SHARED_INFORMATION_MF{1.857047684913544};
  const static constexpr double TEST_MAX_SHARED_INFORMATION_MF{10.973597480731526};
  const static constexpr double TEST_CROSS_ONTOLOGY_SHARED_INFORMATION{0.0};

};

namespace kol = kellerberrin::ontology;

using TestSimilarityClass = kol::TestSharedSimilarity<kol::InformationAncestorMean, AncestorMeanSharedValues>;

BOOST_FIXTURE_TEST_SUITE(TestAncestorMeanSharedSuite, TestSimilarityClass)

#include "kol_test_sharedinfo.h"

BOOST_AUTO_TEST_SUITE_END()
