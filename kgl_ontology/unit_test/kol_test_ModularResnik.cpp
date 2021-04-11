//
// Created by kellerberrin on 5/4/21.
//

#include "kol_test.h"

namespace kol = kellerberrin::ontology;

using ModularTestClass = kol::TestModTermSimilarity<kol::ResnikSimilarity, ModularResnik>;

BOOST_FIXTURE_TEST_SUITE(TestModularResnickSuite, ModularTestClass)

#include "kol_test_modularsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
