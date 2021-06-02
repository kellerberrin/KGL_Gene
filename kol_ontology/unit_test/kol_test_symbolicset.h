//
// Created by kellerberrin on 2/6/21.
//

#ifndef KOL_TEST_SYMBOLICSET_H
#define KOL_TEST_SYMBOLICSET_H
///////////////////////////////////////////////////////////////
// Include this file between the test suite markers
//
// BOOST_FIXTURE_TEST_SUITE(TestAllPairsMaxSetSuite, TestSimilarityClass)
// #include "ggt_test_setsimilarity.h"
// BOOST_AUTO_TEST_SUITE_END()
//
///////////////////////////////////////////////////////////////


#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>


/////////////////////////////////////////
// Basic test
/////////////////////////////////////////


namespace kol = kellerberrin::ontology;

BOOST_AUTO_TEST_CASE(test_set_similarity_empty_sets)
{

  kol::OntologySetType<std::string> termsA;
  kol::OntologySetType<std::string> termsB;
  if (setSimilarity().calculateSimilarity(termsA, termsB) != getValue().TEST_SET_SIMILARITY_EMPTY_SETS ) {

    BOOST_FAIL("Non-zero value on empty sets");

  }
  BOOST_TEST_MESSAGE( "test_set_similarity_empty_sets ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_1_empty_1_good)
{
  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("A0A0B4J269");
  kol::OntologySetType<std::string> termsB;
  if (setSimilarity().calculateSimilarity(termsA, termsB) != getValue().TEST_SET_SIMILARITY_1_EMPTY_1_GOOD) {

    BOOST_FAIL("Non-zero value on empty set");

  }
  BOOST_TEST_MESSAGE( "test_set_similarity_1_empty_1_good ... OK" );

}

///////////////////////////////////////////////////
// Biological process
//////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_BP)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("TUBB3");
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_BP)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("TUBB3");
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_BP)
{

  auto termsA = annotation().getGoTermsForGeneBP("TUBB3");
  auto termsB = annotation().getGoTermsForGeneBP("ARHGAP10");
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_BP ... OK" );

}

///////////////////////////////////////////
// Molecular function
///////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_MF)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneMF("GOPC");
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_MF)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneMF("GOPC");
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_MF)
{

  auto termsA = annotation().getGoTermsForGeneMF("GOPC");
  auto termsB = annotation().getGoTermsForGeneMF("MYO1C");
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_MF ... OK" );

}

/////////////////////////////////////////////////////////////
// Cellular component
/////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_CC)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneCC("TMEM129");
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_CC ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_CC)
{

  kol::OntologySetType<std::string> termsA = annotation().getGoTermsForGeneCC("TMEM129");
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_CC)
{

  auto termsA = annotation().getGoTermsForGeneCC("TMEM129");
  auto termsB = annotation().getGoTermsForGeneCC("TMEM120B");
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_CC ... OK" );

}



#endif //KOL_TEST_SYMBOLICSET_H
