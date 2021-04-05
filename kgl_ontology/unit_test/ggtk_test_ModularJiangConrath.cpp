//
// Created by kellerberrin on 5/4/21.
//



#include "ggtk_test.h"


using ModularTestClass = TestModTermSimilarity<JiangConrathSimilarity, ModularJiangConrath>;

BOOST_FIXTURE_TEST_SUITE(TestModularJiangConrathSuite, ModularTestClass)

#include "ggtk_test_modularsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()

