//
// Created by kellerberrin on 2/4/21.
//
#include "kol_library.h"
#include "kol_test.h"
#include <boost/test/unit_test.hpp>

namespace kellerberrin::ontology {


class TestAnnotateParsers {

public:

  TestAnnotateParsers() = default;

  ~TestAnnotateParsers() = default;


  [[nodiscard]] std::shared_ptr<const TermAnnotation> parseGaf(const PolicyEvidence &policy = PolicyEvidence()) {

    return parseAnnotation(UnitTestDefinitions::gafFileName(), policy);

  }

  [[nodiscard]] std::shared_ptr<const TermAnnotation> parseGene(const PolicyEvidence &policy = PolicyEvidence()) {

    return parseAnnotation(UnitTestDefinitions::geneFileName(), policy);

  }


private:

  [[nodiscard]] std::shared_ptr<const TermAnnotation> parseAnnotation( const std::string &annotation_file,
                                                                       const PolicyEvidence &policy) {

    return ParserAnnotationGaf::parseAnnotationFile(policy, annotation_file);

  }



};

} // namespace

namespace kol = kellerberrin::ontology;

BOOST_FIXTURE_TEST_SUITE(TestAnnotateParsersSuite, kol::TestAnnotateParsers)

//////////////////////////////////////////////////////////
// Parse an annotation file normally with default settings
//////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf)
{

  auto annotation_ptr = parseGaf();
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() == 0 or annotation_ptr->getNumGoTerms() == 0) BOOST_FAIL( "Gaf annotation is empty." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf ... OK" );

}

BOOST_AUTO_TEST_CASE(test_annotation_parser_mgi)
{

  auto annotation_ptr = parseGene();
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() == 0 or annotation_ptr->getNumGoTerms() == 0) BOOST_FAIL( "Gene annotation is empty." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_mgi ... OK" );

}





//////////////////////////////////////////////////
// Bad evidence codes
//////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_bad_custom_evidence_set)
{

  std::vector<kol::GO::EvidenceCode> empty_vector;
  kol::PolicyEvidence bad_policy(empty_vector);
  auto annotation_ptr = parseGaf(bad_policy);
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() != 0 or annotation_ptr->getNumGoTerms() != 0) BOOST_FAIL( "Gaf annotation is non-empty with bad policy." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_bad_custom_evidence_set ... OK" );

}


/////////////////////////////////////////
// Experimental evidence codes
/////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_experimental_evidence_set)
{

  auto experimental_codes = kol::GO::getEvidenceType(kol::GO::EvidenceType::EXPERIMENTAL);
  auto experimental_policy = kol::PolicyEvidence(experimental_codes);
  auto annotation_ptr = parseGaf(experimental_policy);
  BOOST_REQUIRE(annotation_ptr);
  auto gaf_record_vector = annotation_ptr->getAllGAFRecords();
  bool found_bad_codes{false};
  for (auto const& gaf_record_ptr : gaf_record_vector) {

    if (not experimental_policy.isAllowed(gaf_record_ptr->evidenceCode())) {

      found_bad_codes = true;
      break;

    }

  }

  if (found_bad_codes) BOOST_FAIL( "Found an annotated gene with no experimental codes." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_experimental_evidence_set ... OK" );

}




BOOST_AUTO_TEST_SUITE_END()
