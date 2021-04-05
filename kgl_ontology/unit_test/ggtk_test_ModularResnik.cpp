//
// Created by kellerberrin on 5/4/21.
//

#include "ggtk_test.h"

using ModularTestClass = TestModTermSimilarity<ResnikSimilarity, ModularResnik>;

BOOST_FIXTURE_TEST_SUITE(TestModularResnickSuite, ModularTestClass)

#include "ggtk_test_modularsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()
