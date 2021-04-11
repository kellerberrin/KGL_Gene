//
// Created by kellerberrin on 1/4/21.
//


#include <kol_library.h>
#include "kol_test.h"
#include <boost/test/unit_test.hpp>


namespace kellerberrin::ontology {

class TestGoParsers {

public:

  TestGoParsers() = default;

  ~TestGoParsers() = default;

  [[nodiscard]] auto checkOboParser() const {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(go_obo_);

  }

  [[nodiscard]] auto checkXmlParser() const {

    auto xml_parser_ptr = GoParserFactory::createGoParser(GoParserType::XML_GO_STANDARD);
    BOOST_REQUIRE(xml_parser_ptr);
    return xml_parser_ptr->parseGoFile(go_xml_);

  }


  [[nodiscard]] auto checkRelationshipOboParser(const RelationshipPolicyInterface &policy) const {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_ALLOWED, policy);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(go_obo_);

  }

  [[nodiscard]] auto checkRelationshipXmlParser(const RelationshipPolicyInterface &policy) const {

    auto xml_parser_ptr = GoParserFactory::createGoParser(GoParserType::XML_GO_ALLOWED, policy);
    BOOST_REQUIRE(xml_parser_ptr);
    return xml_parser_ptr->parseGoFile(go_xml_);

  }

  [[nodiscard]] bool checkOboFile(const std::string &file_name) const {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->isFileGood(file_name);

  }

  [[nodiscard]] bool checkXmlFile(const std::string &file_name) const {

    auto xml_parser_ptr = GoParserFactory::createGoParser(GoParserType::XML_GO_STANDARD);
    BOOST_REQUIRE(xml_parser_ptr);
    return xml_parser_ptr->isFileGood(file_name);

  }

  [[nodiscard]] const std::string &oboFileName() const { return go_obo_; }

  [[nodiscard]] const std::string &xmlFileName() const { return go_xml_; }

private:

  const std::string go_obo_{UnitTestDefinitions::oboFileName()};
  const std::string go_xml_{UnitTestDefinitions::xmlFileName()};

};

} // namespace

namespace kol = kellerberrin::ontology;

BOOST_FIXTURE_TEST_SUITE(TestGoParserSuite, kol::TestGoParsers)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse go files normally with the standard relationship set
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_parse_obo)
{

  auto graph_ptr = checkOboParser();
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Obo graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_obo ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml)
{

  auto graph_ptr = checkXmlParser();
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Xml graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_xml ... OK" );

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parse a go file with a custom relationship set
//////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_parse_obo_custom_relationships)
{

  auto graph_ptr = checkRelationshipOboParser(kol::StandardRelationshipPolicy());
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Obo graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_obo_custom_relationship ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml_custom_relationships)
{

  auto graph_ptr = checkRelationshipXmlParser(kol::StandardRelationshipPolicy());
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Xml graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_xml_custom_relationship ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_obo_custom_relationships_bad_set)
{

  kol::AllowedSetRelationshipPolicy bad_policy;
  bad_policy.addRelationship("part_of");
  auto graph_ptr = checkRelationshipOboParser(bad_policy);
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() != 0 or graph_ptr->getNumEdges() != 0) BOOST_FAIL( "Obo graph is non-empty." );
  BOOST_TEST_MESSAGE( "test_parse_obo_custom_relationships_bad_set ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml_custom_relationships_bad_set)
{

  kol::AllowedSetRelationshipPolicy bad_policy;
  bad_policy.addRelationship(kol::GO::Relationship::PART_OF);
  auto graph_ptr = checkRelationshipXmlParser(bad_policy);
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() != 0 or graph_ptr->getNumEdges() != 0) BOOST_FAIL( "Xml graph is non-empty." );
  BOOST_TEST_MESSAGE( "test_parse_xml_custom_relationships_bad_set ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_obo_all_relationships)
{

  kol::AllowedSetRelationshipPolicy all_policy;
  all_policy.addRelationship(kol::GO::Relationship::IS_A);
  all_policy.addRelationship(kol::GO::Relationship::PART_OF);
  all_policy.addRelationship(kol::GO::Relationship::REGULATES);
  all_policy.addRelationship(kol::GO::Relationship::NEGATIVELY_REGULATES);
  all_policy.addRelationship(kol::GO::Relationship::POSITIVELY_REGULATES);
  auto graph_ptr = checkRelationshipOboParser(all_policy);
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Obo graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_obo_all_relationships ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml_all_relationships)
{

  kol::AllowedSetRelationshipPolicy all_policy;
  all_policy.addRelationship("is_a");
  all_policy.addRelationship("part_of");
  all_policy.addRelationship("regulates");
  all_policy.addRelationship("positively_regulates");
  all_policy.addRelationship("negatively_regulates");
  auto graph_ptr = checkRelationshipXmlParser(all_policy);
  BOOST_REQUIRE(graph_ptr);
  if( graph_ptr->getNumVertices() == 0 or graph_ptr->getNumEdges() == 0) BOOST_FAIL( "Xml graph is empty." );
  BOOST_TEST_MESSAGE( "test_parse_xml_all_relationships ... OK" );

}

//////////////////////////////////////////////////////////////////////////////
// Nonexistant file and bad format
//////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_parse_obo_bad_file_name)
{

  if( checkOboFile("")) BOOST_FAIL( "Obo file is valid" );
  BOOST_TEST_MESSAGE( "test_parse_obo_bad_file_name ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml_bad_file_name)
{

  if( checkXmlFile("")) BOOST_FAIL( "Xml file is valid" );
  BOOST_TEST_MESSAGE( "test_parse_xml_bad_file_name ... OK" );

}

// Bad Formatting

BOOST_AUTO_TEST_CASE(test_parse_obo_bad_format)
{

  if( checkOboFile(xmlFileName())) BOOST_FAIL( "Obo file is valid" );
  BOOST_TEST_MESSAGE( "test_parse_obo_bad_format ... OK" );

}

BOOST_AUTO_TEST_CASE(test_parse_xml_bad_format)
{

  if( checkXmlFile(oboFileName())) BOOST_FAIL( "Xml file is valid" );
  BOOST_TEST_MESSAGE( "test_parse_xml_bad_format ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
