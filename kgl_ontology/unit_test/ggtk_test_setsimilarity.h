//
// Created by kellerberrin on 5/4/21.
//

#ifndef GGTK_TEST_SETSIMILARITY_H
#define GGTK_TEST_SETSIMILARITY_H
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


BOOST_AUTO_TEST_CASE(test_set_similarity_empty_sets)
{

  OntologySetType<std::string> termsA;
  OntologySetType<std::string> termsB;
  if (setSimilarity().calculateSimilarity(termsA, termsB) != getValue().TEST_SET_SIMILARITY_EMPTY_SETS ) {

    BOOST_FAIL("Non-zero value on empty sets");

  }
  BOOST_TEST_MESSAGE( "test_set_similarity_empty_sets ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_1_empty_1_good)
{
  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("A0A0B4J269", goGraph());
  OntologySetType<std::string> termsB;
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

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("A0A0B4J269", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_BP)
{

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneBP("A0A0B4J269", goGraph());
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_BP)
{

  auto termsA = annotation().getGoTermsForGeneBP("A0A0B4J269", goGraph());
  auto termsB = annotation().getGoTermsForGeneBP("A1A4S6", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_BP ... OK" );

}

///////////////////////////////////////////
// Molecular function
///////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_MF)
{

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneMF("A0A0J9YVX5", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_MF)
{

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneMF("A0A0J9YVX5", goGraph());
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_MF)
{

  auto termsA = annotation().getGoTermsForGeneMF("A0A0J9YVX5", goGraph());
  auto termsB = annotation().getGoTermsForGeneMF("O00159", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_MF ... OK" );

}

/////////////////////////////////////////////////////////////
// Cellular component
/////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_CC)
{

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneCC("A0AVI4", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsA);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_CC ... OK" );

}


BOOST_AUTO_TEST_CASE(test_set_similarity_reflexive_slice_CC)
{

  OntologySetType<std::string> termsA = annotation().getGoTermsForGeneCC("A0AVI4", goGraph());
  auto termsB = termsA;
  if (not termsB.empty()) termsB.erase(termsB.begin());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_REFLEXIVE_SLICE_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_reflexive_slice_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_set_similarity_CC)
{

  auto termsA = annotation().getGoTermsForGeneCC("A0AVI4", goGraph());
  auto termsB = annotation().getGoTermsForGeneCC("A0PK00", goGraph());
  double value = setSimilarity().calculateSimilarity(termsA, termsB);
  BOOST_CHECK_CLOSE(value, getValue().TEST_SET_SIMILARITY_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_set_similarity_CC ... OK" );

}


#endif //KGL_GGTK_TEST_SETSIMILARITY_H
