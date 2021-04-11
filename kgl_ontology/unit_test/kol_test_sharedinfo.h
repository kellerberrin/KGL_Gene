//
// Created by kellerberrin on 6/4/21.
//
///////////////////////////////////////////////////////////////
// Include this file between the test suite markers
//
// BOOST_FIXTURE_TEST_SUITE(TestAdjCoutoSharedSuite, TestSimilarityClass)
// #include "ggt_test_sharedinfo.h"
// BOOST_AUTO_TEST_SUITE_END()
//
///////////////////////////////////////////////////////////////

#ifndef GGTK_TEST_SHAREDINFO_H
#define GGTK_TEST_SHAREDINFO_H

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>



//////////////////////////////////////////////////////////////
// Non-exsitent terms used as input
//////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_shared_information_bad_ids)
{

  if (sharedAnalysis().sharedInformation("bad_id","bad_id2") != getValue().TEST_SHARED_INFORMATION_BAD_IDS ) {

    BOOST_FAIL("Non-zero value on bad id");

  }
  BOOST_TEST_MESSAGE( "test_shared_information_bad_ids ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_1_bad_1_good_id)
{

  if (sharedAnalysis().sharedInformation("GO:0032991","bad_id2") != getValue().TEST_SHARED_INFORMATION_1_BAD_1_GOOD_ID ) {

    BOOST_FAIL("Non-zero value on bad id");

  }
  BOOST_TEST_MESSAGE( "test_shared_information_1_bad_1_good_id ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Deep term
///////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_shared_information_deep_terms_BP)
{

  double value = sharedAnalysis().sharedInformation("GO:0000413", "GO:0043966");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_DEEP_TERMS_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_deep_terms_BP ... OK" );

}

/////////////////////////////////////////////////////////////////////////////////
// Shared Information on CC terms
/////////////////////////////////////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_shared_information_CC_reflexive_sim)
{

  double value = sharedAnalysis().sharedInformation("GO:0043234", "GO:0043234");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_CC_REFLEXIVE_SIM, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_CC_reflexive_sim ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_CC)
{

  double value = sharedAnalysis().sharedInformation("GO:0043234", "GO:0000791");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_CC ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_CC_1_good_1_root)
{

  double value = sharedAnalysis().sharedInformation("GO:0043234", "GO:0005575");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_CC_1_GOOD_1_ROOT, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_CC_1_good_1_root ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_single_term_CC)
{

  double value = sharedAnalysis().sharedInformation("GO:0043234");
  BOOST_CHECK_EQUAL(value, termInformation().getValue("GO:0043234"));
  BOOST_TEST_MESSAGE( "test_shared_information_single_term_CC ... OK" );

}


BOOST_AUTO_TEST_CASE(test_max_shared_information_CC)
{

  double value = sharedAnalysis().maxInformationContent("GO:0043234");
  BOOST_CHECK_CLOSE( value, getValue().TEST_MAX_SHARED_INFORMATION_CC, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_max_shared_information_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////////
// Shared Information in BP
/////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_shared_information_BP)
{

  double value = sharedAnalysis().sharedInformation("GO:0007155", "GO:0044406");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_single_term_BP)
{

  double value = sharedAnalysis().sharedInformation("GO:0007155");
  BOOST_CHECK_EQUAL(value, termInformation().getValue("GO:0007155"));
  BOOST_TEST_MESSAGE( "test_shared_information_single_term_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_max_shared_information_BP)
{

  double value = sharedAnalysis().maxInformationContent("GO:0007155");
  BOOST_CHECK_CLOSE( value, getValue().TEST_MAX_SHARED_INFORMATION_BP, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_max_shared_information_BP ... OK" );

}

/////////////////////////////////////////////////////////////////////////
// Shared Information in MF
/////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_shared_information_MF)
{

  double value = sharedAnalysis().sharedInformation("GO:0051192", "GO:0050662");
  BOOST_CHECK_CLOSE( value, getValue().TEST_SHARED_INFORMATION_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_shared_information_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_shared_information_single_term_MF)
{

  double value = sharedAnalysis().sharedInformation("GO:0051192");
  BOOST_CHECK_EQUAL(value, termInformation().getValue("GO:0051192"));
  BOOST_TEST_MESSAGE( "test_shared_information_single_term_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_max_shared_information_MF)
{

  double value = sharedAnalysis().maxInformationContent("GO:0051192");
  BOOST_CHECK_CLOSE( value, getValue().TEST_MAX_SHARED_INFORMATION_MF, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_max_shared_information_MF ... OK" );

}


//////////////////////////////////////////////////////////////////////////////
// Cross Ontology Shared Information
//////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_cross_ontology_shared_information)
{

  double value = sharedAnalysis().sharedInformation("GO:0051192", "GO:0007155");
  BOOST_CHECK_CLOSE( value, getValue().TEST_CROSS_ONTOLOGY_SHARED_INFORMATION, TEST_ACCURACY_PERCENT);
  BOOST_TEST_MESSAGE( "test_cross_ontology_shared_information ... OK" );

}



#endif //GGTK_TEST_SHAREDINFO_H
