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


  [[nodiscard]] std::unique_ptr<AnnotationData> parseGaf(const PolicyAllowedEvidence &policy = PolicyAllowedEvidence()) {

    return parseAnnotation(AnnotationParserType::GAF_ANNO_PARSER, UnitTestDefinitions::gafFileName(), policy);

  }

  [[nodiscard]] std::unique_ptr<AnnotationData> parseEntrez(const PolicyAllowedEvidence &policy = PolicyAllowedEvidence()) {

    return parseAnnotation(AnnotationParserType::ENTREZ_ANNO_PARSER, UnitTestDefinitions::entrezFileName(), policy);

  }

  [[nodiscard]] std::unique_ptr<AnnotationData> parseGene(const PolicyAllowedEvidence &policy = PolicyAllowedEvidence()) {

    return parseAnnotation(AnnotationParserType::MGI_ANNO_PARSER, UnitTestDefinitions::geneFileName(), policy);

  }

  [[nodiscard]] bool verifyGaf(const std::string &file_name = UnitTestDefinitions::gafFileName()) {

    return verifyAnnotation(AnnotationParserType::GAF_ANNO_PARSER, file_name);

  }

  [[nodiscard]] bool verifyEntrez(const std::string &file_name = UnitTestDefinitions::entrezFileName()) {

    return verifyAnnotation(AnnotationParserType::ENTREZ_ANNO_PARSER, file_name);

  }

  [[nodiscard]] bool verifyGene(const std::string &file_name = UnitTestDefinitions::geneFileName()) {

    return verifyAnnotation(AnnotationParserType::MGI_ANNO_PARSER, file_name);

  }


private:

  [[nodiscard]] std::unique_ptr<AnnotationData> parseAnnotation(AnnotationParserType parser_type,
                                                                const std::string &annotation_file,
                                                                const PolicyAllowedEvidence &policy) {

    auto anno_parser_ptr = ParserAnnotationFactory::createAnnotationParser(parser_type, policy);
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(annotation_file);

  }

  [[nodiscard]] bool verifyAnnotation(AnnotationParserType parser_type,
                                      const std::string &annotation_file) {

    auto anno_parser_ptr = ParserAnnotationFactory::createAnnotationParser(parser_type, PolicyAllowedEvidence());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->isFileGood(annotation_file);

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

BOOST_AUTO_TEST_CASE(test_annotation_parser_entrez)
{

  auto annotation_ptr = parseEntrez();
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() == 0 or annotation_ptr->getNumGoTerms() == 0) BOOST_FAIL( "Entrez annotation is empty." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_entrez ... OK" );

}

//////////////////////////////////////////////////////
// Test Parser with bad input
//////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_nonexistent_file)
{

  if ( verifyGaf("") ) BOOST_FAIL( "Empty Gaf file is verified valid" );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_nonexistent_file ... OK" );

}


BOOST_AUTO_TEST_CASE(test_annotation_parser_mgi_nonexistent_file)
{

  if ( verifyGene("") ) BOOST_FAIL( "Empty Gene Annotation file is verified valid" );
  BOOST_TEST_MESSAGE( "test_annotation_parser_mgi_nonexistent_file ... OK" );

}



BOOST_AUTO_TEST_CASE(test_annotation_parser_entrez_nonexistent_file)
{

  if ( verifyEntrez("") ) BOOST_FAIL( "Empty Entrez file is verified valid" );
  BOOST_TEST_MESSAGE( "test_annotation_parser_entrez_nonexistent_file ... OK" );

}

//////////////////////////////////////////////////////
// Bad file format
//////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_bad_format)
{

  if ( verifyGaf(UnitTestDefinitions::entrezFileName()) ) BOOST_FAIL( "Invalid Gaf file is verified valid" );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_bad_format ... OK" );

}


BOOST_AUTO_TEST_CASE(test_annotation_parser_entrez_bad_format)
{

  if ( verifyEntrez(UnitTestDefinitions::gafFileName()) ) BOOST_FAIL( "Invalid Entrez Annotation file is verified valid" );
  BOOST_TEST_MESSAGE( "test_annotation_parser_entrez_bad_format ... OK" );

}

//////////////////////////////////////////////////
// Bad evidence codes
//////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_bad_custom_evidence_set)
{

  std::vector<kol::GO::EvidenceCode> empty_vector;
  kol::PolicyAllowedEvidence bad_policy(empty_vector);
  auto annotation_ptr = parseGaf(bad_policy);
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() != 0 or annotation_ptr->getNumGoTerms() != 0) BOOST_FAIL( "Gaf annotation is non-empty with bad policy." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_bad_custom_evidence_set ... OK" );

}

BOOST_AUTO_TEST_CASE(test_annotation_parser_entrez_bad_custom_evidence_set)
{

  std::vector<kol::GO::EvidenceCode> empty_vector;
  kol::PolicyAllowedEvidence bad_policy(empty_vector);
  auto annotation_ptr = parseEntrez(bad_policy);
  BOOST_REQUIRE(annotation_ptr);
  if ( annotation_ptr->getNumGenes() != 0 or annotation_ptr->getNumGoTerms() != 0) BOOST_FAIL( "Entrez annotation is non-empty with bad policy." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_entrez_bad_custom_evidence_set ... OK" );

}

/////////////////////////////////////////
// Experimental evidence codes
/////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_annotation_parser_gaf_experimental_evidence_set)
{

  auto experimental_codes = kol::GO::getEvidenceType(kol::GO::EvidenceType::EXPERIMENTAL);
  auto experimental_policy = kol::PolicyAllowedEvidence(experimental_codes);
  auto annotation_ptr = parseGaf(experimental_policy);
  BOOST_REQUIRE(annotation_ptr);
  auto gene_list = annotation_ptr->getAllGenes();

  bool no_experimental_codes{false};
  for (auto const& gene : gene_list) {

    auto code_list = annotation_ptr->getGoTermsEvidenceForGene(gene);
    bool found_experimental_codes{false};
    for (auto const& code_text : code_list) {

      if (experimental_policy.isAllowed(kol::GO::evidenceStringToCode(code_text))) {

        found_experimental_codes = true;
        break;

      }

    }
    if (not found_experimental_codes) {

      no_experimental_codes = true;
      break;

    }

  }
  if (no_experimental_codes) BOOST_FAIL( "Found an annotated gene with no experimental codes." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_gaf_experimental_evidence_set ... OK" );

}


BOOST_AUTO_TEST_CASE(test_annotation_parser_entrez_experimental_evidence_set)
{


  auto experimental_codes = kol::GO::getEvidenceType(kol::GO::EvidenceType::EXPERIMENTAL);
  auto experimental_policy = kol::PolicyAllowedEvidence(experimental_codes);
  auto annotation_ptr = parseEntrez(experimental_policy);
  BOOST_REQUIRE(annotation_ptr);
  auto gene_list = annotation_ptr->getAllGenes();

  bool no_experimental_codes{false};
  for (auto const& gene : gene_list) {

    auto code_list = annotation_ptr->getGoTermsEvidenceForGene(gene);
    bool found_experimental_codes{false};
    for (auto const& code_text : code_list) {

      if (experimental_policy.isAllowed(kol::GO::evidenceStringToCode(code_text))) {

        found_experimental_codes = true;
        break;

      }

    }
    if (not found_experimental_codes) {

      no_experimental_codes = true;
      break;

    }

  }
  if (no_experimental_codes) BOOST_FAIL( "Found an annotated gene with no experimental codes." );
  BOOST_TEST_MESSAGE( "test_annotation_parser_entrez_experimental_evidence_set ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
