//
// Created by kellerberrin on 19/5/21.
//


#include "kol_test.h"


namespace kol = kellerberrin::ontology;

using ModularTestClass = kol::TestModTermSimilarity<kol::RelevanceSimilarity, kol::ModularRelevance>;

BOOST_FIXTURE_TEST_SUITE(TestModularRelevanceSuite, ModularTestClass)

#include "kol_test_modularsimilarity.h"

BOOST_AUTO_TEST_SUITE_END()

