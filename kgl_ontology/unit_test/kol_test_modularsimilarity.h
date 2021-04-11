//
// Created by kellerberrin on 5/4/21.
//
///////////////////////////////////////////////////////////////
// Include this file between the test suite markers
//
// BOOST_FIXTURE_TEST_SUITE(TestModularJiangConrathSuite, ModularTestClass)
// #include "ggt_test_modsimilarity.h"
// BOOST_AUTO_TEST_SUITE_END()
//
///////////////////////////////////////////////////////////////


#ifndef GGTK_TEST_MODSIMILARITY_H
#define GGTK_TEST_MODSIMILARITY_H

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_bad_ids)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("bad_id","bad_id2");
  if (termSimilarity().calculateTermSimilarity("bad_id","bad_id2") !=  modular_similarity) BOOST_FAIL("Non-zero value on bad id");
  BOOST_TEST_MESSAGE( "test_similarity_bad_ids ... OK" );

}


BOOST_AUTO_TEST_CASE(test_similarity_1_bad_1_good_id)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("bad_id","bad_id2");
  if (termSimilarity().calculateTermSimilarity("GO:0032991","bad_id2") != modular_similarity ) BOOST_FAIL("Non-zero value on bad id");
  BOOST_TEST_MESSAGE( "test_similarity_1_bad_1_good_id ... OK" );

}


////////////////////////////////////////////////////////
// Similarity on CC terms
///////////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_similarity_CC_reflexive_sim)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0043234", "GO:0043234");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0043234", "GO:0043234");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_similarity_CC_reflexive_sim ... OK" );

}


BOOST_AUTO_TEST_CASE(test_similarity_CC)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0043234", "GO:0000791");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0043234", "GO:0000791");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_similarity_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_similarity_CC_1_good_1_root)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0043234", "GO:0005575");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0043234", "GO:0005575");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_similarity_CC_1_good_1_root ... OK" );

}


/////////////////////////////////////////////////////////////////////////
// Normalized Similarity [0-1], on CC terms
/////////////////////////////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_normalized_similarity_CC_reflexive_sim)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0043234");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0043234");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_CC_reflexive_sim ... OK" );

}


BOOST_AUTO_TEST_CASE(test_normalized_similarity_CC)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0000791");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0000791");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_CC_1_good_1_root)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0005575");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0043234", "GO:0005575");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_CC_1_good_1_root ... OK" );

}


///////////////////////////////////////////////////////////////////////
// Similarity in BP
///////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_similarity_BP)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0007155", "GO:0044406");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0007155", "GO:0044406");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_similarity_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_BP)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0007155", "GO:0044406");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0007155", "GO:0044406");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_BP ... OK" );

}


//////////////////////////////////////////////////////////////////////////
// Similarity in MF
/////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_similarity_MF)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0051192", "GO:0050662");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0051192", "GO:0050662");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_similarity_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_normalized_similarity_MF)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0051192", "GO:0050662");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0051192", "GO:0050662");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_normalized_similarity_MF ... OK" );

}

////////////////////////////////////////////////////////////////////////////////
// Cross Ontology similarity
////////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_cross_ontology_similarity)
{

  double modular_similarity = termModSimilarity().calculateTermSimilarity("GO:0051192", "GO:0007155");
  double similarity = termSimilarity().calculateTermSimilarity("GO:0051192", "GO:0007155");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_cross_ontology_similarity ... OK" );

}

BOOST_AUTO_TEST_CASE(test_cross_ontology_normalized_similarity)
{

  double modular_similarity = termModSimilarity().calculateNormalizedTermSimilarity("GO:0051192", "GO:0007155");
  double similarity = termSimilarity().calculateNormalizedTermSimilarity("GO:0051192", "GO:0007155");
  BOOST_CHECK_CLOSE( similarity, modular_similarity, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_cross_ontology_normalized_similarity ... OK" );

}


#endif //GGTK_TEST_MODSIMILARITY_H
