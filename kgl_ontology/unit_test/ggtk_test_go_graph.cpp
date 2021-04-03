//
// Created by kellerberrin on 1/4/21.
//
#include <ggtk.h>
#include "ggtk_test.h"
#include <boost/test/unit_test.hpp>

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestGoGraph
{

public:

  TestGoGraph() = default;
  ~TestGoGraph() = default;

  [[nodiscard]] static const GoGraph& goGraph() {

    if (not static_graph_) {

      static_graph_ = getGoGraph();

    }
    BOOST_REQUIRE(static_graph_);
    return *static_graph_;

  }

  // Utility function.
  // Convert an OntologySetType (which may be a std::unordered_set<>) into a std::set<> for convenient '==' comparison.
  template<class T>
  [[nodiscard]] static std::set<T> convertSet(OntologySetType<T>&& from_set) {

    if constexpr(std::is_same<std::set<T>, OntologySetType<T>>::value) {

      return std::set<T>(from_set);

    } else {

      std::set<T> plain_set;
      for (auto&& element : from_set) {

        plain_set.insert(element);

      }
      return plain_set;

    }

  }

private:

  [[nodiscard]] static std::unique_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  inline static std::unique_ptr<const GoGraph> static_graph_;

};


BOOST_FIXTURE_TEST_SUITE(TestGoGraphSuite, TestGoGraph)

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertex and Edge count accessors
//////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_vertex_count_accessor)
{

  const size_t vertex_count{42718};
  if (goGraph().getNumVertices() != vertex_count) BOOST_FAIL("Obo graph vertex count is incorrect" );
  BOOST_TEST_MESSAGE( "test_vertex_count_accessor ... OK" );

}

BOOST_AUTO_TEST_CASE(test_edge_count_accessor)
{

  const size_t edge_count{81084};
  if (goGraph().getNumEdges() != edge_count) BOOST_FAIL("Obo graph edge count is incorrect" );
  BOOST_TEST_MESSAGE( "test_edge_count_accessor ... OK" );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Term attribute access, Name Description and Ontology type
//////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_name_bad_id)
{

  const std::string bad_term("GO:00507");
  if (not goGraph().getTermName(bad_term).empty()) BOOST_FAIL( "Bad term found" );
  BOOST_TEST_MESSAGE( "test_term_name_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_bad_id)
{

  const std::string bad_term("GO:00507");
  if (not goGraph().getTermDescription(bad_term).empty()) BOOST_FAIL( "Bad term description found" );
  BOOST_TEST_MESSAGE( "test_term_desc_bad_id ... OK" );

}

///////////////////////////////////////////////////////
// Vertex and Edge count accessors
//////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_name_BP)
{

  const std::string BP_term("GO:0050789");
  const std::string BP_name("regulation of biological process");
  if (goGraph().getTermName(BP_term) != BP_name) BOOST_FAIL( "BP term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_BP)
{

  const std::string BP_term("GO:0050789");
  const std::string part_desc("Any process that modulates the frequency, rate or extent of a biological process");
  if (goGraph().getTermDescription(BP_term).find(part_desc) == std::string::npos) BOOST_FAIL( "BP description for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_name_MF)
{

  const std::string MF_term("GO:0005385");
  const std::string MF_name("zinc ion transmembrane transporter activity");
  if (goGraph().getTermName(MF_term) != MF_name) BOOST_FAIL( "MF term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_MF)
{

  const std::string MF_term("GO:0005385");
  const std::string part_desc("\"Enables the transfer of zinc (Zn) ions from one side of a membrane to the other.\" [GOC:dgf]");
  if (goGraph().getTermDescription(MF_term).find(part_desc) == std::string::npos) BOOST_FAIL( "MF description for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_name_CC)
{

  const std::string CC_term("GO:0009898");
  const std::string CC_name("cytoplasmic side of plasma membrane");
  if (goGraph().getTermName(CC_term) != CC_name) BOOST_FAIL( "CC term for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_CC)
{

  const std::string CC_term("GO:0009898");
  const std::string part_desc("\"The leaflet the plasma membrane that faces the cytoplasm and any proteins embedded or anchored in it or attached to its surface.\" [GOC:dos, GOC:tb]");
  if (goGraph().getTermDescription(CC_term).find(part_desc) == std::string::npos) BOOST_FAIL( "CC description for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Ontology type and code
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_ontology_bad_id)
{

  const std::string bad_term("GO:00098");
  if (GO::ontologyToString(goGraph().getTermOntology(bad_term)) != GO::ONTOLOGY_ERROR_TEXT) BOOST_FAIL( "Bad ontology term for: " + bad_term + " found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_bad_id)
{

  const std::string bad_term("GO:00098");
  if (goGraph().getTermOntology(bad_term) != GO::Ontology::ONTO_ERROR) BOOST_FAIL( "Bad ontology term for: " + bad_term + " found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_bad_id ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_BP)
{

  const std::string BP_term("GO:0022403");
  if (GO::ontologyToString(goGraph().getTermOntology(BP_term)) != GO::ONTOLOGY_BIOLOGICAL_PROCESS_TEXT) BOOST_FAIL( "BP ontology term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_BP)
{

  const std::string BP_term("GO:0022403");
  if (goGraph().getTermOntology(BP_term) != GO::Ontology::BIOLOGICAL_PROCESS) BOOST_FAIL( "BP ontology term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_MF)
{

  const std::string MF_term("GO:0048037");
  if (GO::ontologyToString(goGraph().getTermOntology(MF_term)) != GO::ONTOLOGY_MOLECULAR_FUNCTION_TEXT) BOOST_FAIL( "MF ontology term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_MF)
{

  const std::string MF_term("GO:0048037");
  if (goGraph().getTermOntology(MF_term) != GO::Ontology::MOLECULAR_FUNCTION) BOOST_FAIL( "MF ontology term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_CC)
{

  const std::string CC_term("GO:0005911");
  if (GO::ontologyToString(goGraph().getTermOntology(CC_term)) != GO::ONTOLOGY_CELLULAR_COMPONENT_TEXT) BOOST_FAIL( "CC ontology term for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_CC)
{

  const std::string CC_term("GO:0005911");
  if (goGraph().getTermOntology(CC_term) != GO::Ontology::CELLULAR_COMPONENT) BOOST_FAIL( "CC ontology term for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Relative accessors, Ancestors, Descendants, Parents, Children
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////
// Ancestors
/////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_ancestors_bad_id)
{

  const std::string Bad_term("GO:00098");
  if (not goGraph().getAncestorTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found ancestors" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_BP)
{

  const std::string BP_term("GO:0022403");
  const std::set<std::string> ancestors{"GO:0044848", "GO:0008150"};
  if (convertSet(goGraph().getAncestorTerms(BP_term)) != ancestors) BOOST_FAIL( "BP term: " + BP_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_MF)
{

  const std::string MF_term("GO:0048037");
  const std::set<std::string> ancestors{"GO:0005488", "GO:0003674"};
  if (convertSet(goGraph().getAncestorTerms(MF_term)) != ancestors) BOOST_FAIL( "MF term: " + MF_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_CC)
{

  const std::string CC_term("GO:0005911");
  const std::set<std::string> ancestors{"GO:0005575", "GO:0030054"};
  if (convertSet(goGraph().getAncestorTerms(CC_term)) != ancestors) BOOST_FAIL( "CC term: " + CC_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_CC ... OK" );

}


/////////////////////////////////////////////
// Descendants
/////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_term_descendants_bad_id)
{

  const std::string Bad_term("GO:00098");
  if (not goGraph().getDescendantTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found descendants" );
  BOOST_TEST_MESSAGE( "test_term_descendants_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_descendants_BP)
{

  const std::string BP_term("GO:0051318");
  const std::set<std::string> descendants{"GO:0000080", "GO:0051330"};
  if (convertSet(goGraph().getDescendantTerms(BP_term)) != descendants) BOOST_FAIL( "BP term: " + BP_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_descendants_MF)
{

  const std::string MF_term("GO:0042165");
  const std::set<std::string> descendants{"GO:0031626", "GO:0042166"};
  if (convertSet(goGraph().getDescendantTerms(MF_term)) != descendants) BOOST_FAIL( "MF term: " + MF_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_descendants_CC)
{
  const std::string CC_term("GO:0030057");
  const std::set<std::string> descendants{"GO:0090635", "GO:0090636", "GO:0090637"};
  if (convertSet(goGraph().getDescendantTerms(CC_term)) != descendants) BOOST_FAIL( "CC term: " + CC_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_CC ... OK" );

}

////////////////////////////////////////////////
// Parents, immediate ancestors
////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_term_parents_bad_id)
{

  const std::string Bad_term("GO:0022");
  if (not goGraph().getParentTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found parents" );
  BOOST_TEST_MESSAGE( "test_term_parents_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_BP)
{

  const std::string BP_term("GO:0022403");
  const std::set<std::string> parents{"GO:0044848"};
  if (convertSet(goGraph().getParentTerms(BP_term)) != parents) BOOST_FAIL( "BP term: " + BP_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_MF)
{

  const std::string MF_term("GO:0048037");
  const std::set<std::string> parents{"GO:0005488"};
  if (convertSet(goGraph().getParentTerms(MF_term)) != parents) BOOST_FAIL( "MF term: " + MF_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_CC)
{

  const std::string CC_term("GO:0005911");
  const std::set<std::string> parents{"GO:0030054"};
  if (convertSet(goGraph().getParentTerms(CC_term)) != parents) BOOST_FAIL( "CC term: " + CC_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_CC ... OK" );

}

////////////////////////////////////////////////
// Children, immediate descendants
////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_term_children_bad_id)
{

  const std::string Bad_term("GO:00098");
  if (not goGraph().getChildTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found children" );
  BOOST_TEST_MESSAGE( "test_term_children_bad_id ... OK" );

}
BOOST_AUTO_TEST_CASE(test_term_children_BP)
{

  const std::string BP_term("GO:0051319");
  const std::set<std::string> children{"GO:0051331", "GO:0000085"};
  if (convertSet(goGraph().getChildTerms(BP_term)) != children) BOOST_FAIL( "BP term: " + BP_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_children_MF)
{
  const std::string MF_term("GO:0001071");
  const std::set<std::string> children{"GO:0003700", "GO:0001070"};
  if (convertSet(goGraph().getChildTerms(MF_term)) != children) BOOST_FAIL( "MF term: " + MF_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_children_CC)
{

  const std::string CC_term("GO:0031974");
  const std::set<std::string> children{"GO:1904724", "GO:1904813", "GO:0043233", "GO:0031970"};
  if (convertSet(goGraph().getChildTerms(CC_term)) != children) BOOST_FAIL( "CC term: " + CC_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////
// Retrieve terms sets and ontology term sets
/////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_all_terms_accessor)
{

  const size_t term_count = 42718;
  if (goGraph().getAllTerms().size() != term_count) BOOST_FAIL( "All terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_BP)
{

  const size_t term_count = 28651;
  if (goGraph().getAllTermsBP().size() != term_count) BOOST_FAIL( "All BP terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_MF)
{

  const size_t term_count = 10160;
  if (goGraph().getAllTermsMF().size() != term_count) BOOST_FAIL( "All MF terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_CC)
{

  const size_t term_count = 3907;
  if (goGraph().getAllTermsCC().size() != term_count) BOOST_FAIL( "All CC terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_CC ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Test filtering functions
///////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_filter_for_BP)
{

  OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0033631", "GO:0033627"};
  if (convertSet(goGraph().filterSetForBP(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all BP");
  BOOST_TEST_MESSAGE( "test_filter_for_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_filter_for_MF)
{

  OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0004099", "GO:0019213"};
  if (convertSet(goGraph().filterSetForMF(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all MF");
  BOOST_TEST_MESSAGE( "test_filter_for_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_filter_for_CC)
{

  OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0044225", "GO:0043025"};
  if (convertSet(goGraph().filterSetForCC(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all CC");
  BOOST_TEST_MESSAGE( "test_filter_for_CC ... OK" );

}

///////////////////////////////////////////////////////////
// Test accessors for the root terms of the 3 ontologies
//////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_root_term_access_BP)
{

  if (goGraph().getTermRoot("GO:0033631") != GO::getRootTermBP()) BOOST_FAIL( "Failed to get root term for BP");
  BOOST_TEST_MESSAGE( "test_root_term_access_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_term_access_MF)
{

  if (goGraph().getTermRoot("GO:0019213") != GO::getRootTermMF()) BOOST_FAIL( "Failed to get root term for MF");
  BOOST_TEST_MESSAGE( "test_root_term_access_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_term_access_CC)
{

  if (goGraph().getTermRoot("GO:0044225") != GO::getRootTermCC()) BOOST_FAIL( "Failed to get root term for CC");
  BOOST_TEST_MESSAGE( "test_root_term_access_CC ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
