//
// Created by kellerberrin on 2/6/21.
//

#include "kol_test.h"

struct SymbolicSetValues {

  const static constexpr double TEST_SET_SIMILARITY_EMPTY_SETS{0.0};
  const static constexpr double TEST_SET_SIMILARITY_1_EMPTY_1_GOOD{0.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_BP{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP{0.93243503607289213};
  const static constexpr double TEST_SET_SIMILARITY_BP{0.48639610351129914};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_MF{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF{0.8877085699023004};
  const static constexpr double TEST_SET_SIMILARITY_MF{0.65884469817959812};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_CC{1.0};
  const static constexpr double TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC{0.95167796503881386};
  const static constexpr double TEST_SET_SIMILARITY_CC{0.83620454724454052};


};

namespace kol = kellerberrin::ontology;

using TestSymbolicClass = kol::TestSetTermSymbolicSimilarity<kol::SimilarityJiangConrath, kol::SetSimilarityBestMatchAverage, SymbolicSetValues>;

BOOST_FIXTURE_TEST_SUITE(SymbolicSetSuite, TestSymbolicClass)

#include "kol_test_symbolicset.h"

BOOST_AUTO_TEST_SUITE_END()
