//
// Created by kellerberrin on 1/4/21.
//
#include <kol_library.h>
#include "kol_test.h"
#include "contrib/kol_GoGraphImpl.h"
#include <boost/test/unit_test.hpp>

namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestGoGraph {

public:

  TestGoGraph() = default;

  ~TestGoGraph() = default;

  [[nodiscard]] static const GoGraph &goGraph() {

    if (not static_graph_) {

      static_graph_ = getGoGraph();

    }
//    BOOST_REQUIRE(static_graph_);
    return *static_graph_;

  }

  // Utility function.
  // Convert an OntologySetType (which may be a std::unordered_set<>) into a std::set<> for convenient '==' comparison.
  template<class T>
  [[nodiscard]] static std::set<T> convertSet(OntologySetType<T> &&from_set) {

    if constexpr(std::is_same<std::set<T>, OntologySetType<T>>::value) {

      return std::set<T>(from_set);

    } else {

      std::set<T> plain_set;
      for (auto &&element : from_set) {

        plain_set.insert(element);

      }
      return plain_set;

    }

  }

private:

  [[nodiscard]] static std::shared_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  inline static std::shared_ptr<const GoGraph> static_graph_;

};

} // namespace


namespace kol = kellerberrin::ontology;


BOOST_FIXTURE_TEST_SUITE(TestGoGraphSuite, kol::TestGoGraph)

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vertex and Edge count accessors
//////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_vertex_count_accessor)
{

  const size_t vertex_count{42718};
  if (goGraph().getGoGraphImpl().getNumVertices() != vertex_count) BOOST_FAIL("Obo graph vertex count is incorrect" );
  BOOST_TEST_MESSAGE( "test_vertex_count_accessor ... OK" );

}

BOOST_AUTO_TEST_CASE(test_edge_count_accessor)
{

  const size_t edge_count{81084};
  if (goGraph().getGoGraphImpl().getNumEdges() != edge_count) BOOST_FAIL("Obo graph edge count is incorrect" );
  BOOST_TEST_MESSAGE( "test_edge_count_accessor ... OK" );

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Term attribute access, Name Description and Ontology type
//////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_name_bad_id)
{

  const std::string bad_term("GO:00507");
  if (not goGraph().getGoGraphImpl().getTermName(bad_term).empty()) BOOST_FAIL( "Bad term found" );
  BOOST_TEST_MESSAGE( "test_term_name_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_bad_id)
{

  const std::string bad_term("GO:00507");
  if (not goGraph().getGoGraphImpl().getTermDescription(bad_term).empty()) BOOST_FAIL( "Bad term description found" );
  BOOST_TEST_MESSAGE( "test_term_desc_bad_id ... OK" );

}

///////////////////////////////////////////////////////
// Vertex and Edge count accessors
//////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_name_BP)
{

  const std::string BP_term("GO:0050789");
  const std::string BP_name("regulation of biological process");
  if (goGraph().getGoGraphImpl().getTermName(BP_term) != BP_name) BOOST_FAIL( "BP term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_BP)
{

  const std::string BP_term("GO:0050789");
  const std::string part_desc("Any process that modulates the frequency, rate or extent of a biological process");
  if (goGraph().getGoGraphImpl().getTermDescription(BP_term).find(part_desc) == std::string::npos) BOOST_FAIL( "BP description for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_name_MF)
{

  const std::string MF_term("GO:0005385");
  const std::string MF_name("zinc ion transmembrane transporter activity");
  if (goGraph().getGoGraphImpl().getTermName(MF_term) != MF_name) BOOST_FAIL( "MF term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_MF)
{

  const std::string MF_term("GO:0005385");
  const std::string part_desc("\"Enables the transfer of zinc (Zn) ions from one side of a membrane to the other.\" [GOC:dgf]");
  if (goGraph().getGoGraphImpl().getTermDescription(MF_term).find(part_desc) == std::string::npos) BOOST_FAIL( "MF description for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_name_CC)
{

  const std::string CC_term("GO:0009898");
  const std::string CC_name("cytoplasmic side of plasma membrane");
  if (goGraph().getGoGraphImpl().getTermName(CC_term) != CC_name) BOOST_FAIL( "CC term for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_name_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_desc_CC)
{

  const std::string CC_term("GO:0009898");
  const std::string part_desc("\"The leaflet the plasma membrane that faces the cytoplasm and any proteins embedded or anchored in it or attached to its surface.\" [GOC:dos, GOC:tb]");
  if (goGraph().getGoGraphImpl().getTermDescription(CC_term).find(part_desc) == std::string::npos) BOOST_FAIL( "CC description for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_desc_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Ontology type and code
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_term_ontology_bad_id)
{

  const std::string bad_term("GO:00098");
  if (kol::GO::ontologyToString(goGraph().getGoGraphImpl().getTermOntology(bad_term)) != kol::GO::ONTOLOGY_ERROR_TEXT) BOOST_FAIL( "Bad ontology term for: " + bad_term + " found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_bad_id)
{

  const std::string bad_term("GO:00098");
  if (goGraph().getGoGraphImpl().getTermOntology(bad_term) != kol::GO::Ontology::ONTO_ERROR) BOOST_FAIL( "Bad ontology term for: " + bad_term + " found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_bad_id ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_BP)
{

  const std::string BP_term("GO:0022403");
  if (kol::GO::ontologyToString(goGraph().getGoGraphImpl().getTermOntology(BP_term)) != kol::GO::ONTOLOGY_BIOLOGICAL_PROCESS_TEXT) BOOST_FAIL( "BP ontology term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_BP)
{

  const std::string BP_term("GO:0022403");
  if (goGraph().getGoGraphImpl().getTermOntology(BP_term) != kol::GO::Ontology::BIOLOGICAL_PROCESS) BOOST_FAIL( "BP ontology term for: " + BP_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_MF)
{

  const std::string MF_term("GO:0048037");
  if (kol::GO::ontologyToString(goGraph().getGoGraphImpl().getTermOntology(MF_term)) != kol::GO::ONTOLOGY_MOLECULAR_FUNCTION_TEXT) BOOST_FAIL( "MF ontology term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_MF)
{

  const std::string MF_term("GO:0048037");
  if (goGraph().getGoGraphImpl().getTermOntology(MF_term) != kol::GO::Ontology::MOLECULAR_FUNCTION) BOOST_FAIL( "MF ontology term for: " + MF_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_code_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_ontology_CC)
{

  const std::string CC_term("GO:0005911");
  if (kol::GO::ontologyToString(goGraph().getGoGraphImpl().getTermOntology(CC_term)) != kol::GO::ONTOLOGY_CELLULAR_COMPONENT_TEXT) BOOST_FAIL( "CC ontology term for: " + CC_term + " not found" );
  BOOST_TEST_MESSAGE( "test_term_ontology_CC ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ontology_code_CC)
{

  const std::string CC_term("GO:0005911");
  if (goGraph().getGoGraphImpl().getTermOntology(CC_term) != kol::GO::Ontology::CELLULAR_COMPONENT) BOOST_FAIL( "CC ontology term for: " + CC_term + " not found" );
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
  if (not goGraph().getGoGraphImpl().getAncestorTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found ancestors" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_BP)
{

  const std::string BP_term("GO:0022403");
  const std::set<std::string> ancestors{"GO:0044848", "GO:0008150"};
  if (convertSet(goGraph().getGoGraphImpl().getAncestorTerms(BP_term)) != ancestors) BOOST_FAIL( "BP term: " + BP_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_MF)
{

  const std::string MF_term("GO:0048037");
  const std::set<std::string> ancestors{"GO:0005488", "GO:0003674"};
  if (convertSet(goGraph().getGoGraphImpl().getAncestorTerms(MF_term)) != ancestors) BOOST_FAIL( "MF term: " + MF_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_ancestors_CC)
{

  const std::string CC_term("GO:0005911");
  const std::set<std::string> ancestors{"GO:0005575", "GO:0030054"};
  if (convertSet(goGraph().getGoGraphImpl().getAncestorTerms(CC_term)) != ancestors) BOOST_FAIL( "CC term: " + CC_term + " ancestors not found" );
  BOOST_TEST_MESSAGE( "test_term_ancestors_CC ... OK" );

}


/////////////////////////////////////////////
// Descendants
/////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_term_descendants_bad_id)
{

  const std::string Bad_term("GO:00098");
  if (not goGraph().getGoGraphImpl().getDescendantTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found descendants" );
  BOOST_TEST_MESSAGE( "test_term_descendants_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_terms_with_zero_descendants)
{

  const size_t leaf_terms{23717};
  auto const term_set = goGraph().getGoGraphImpl().getAllTerms();
  size_t leaf_count{0};
  for (auto const& term : term_set) {

    if (goGraph().getGoGraphImpl().isLeaf(term)) {

      ++leaf_count;

    }

  }
  if (leaf_count != leaf_terms) BOOST_FAIL( "Unexpected leaf terms: " + std::to_string(leaf_count) );
  BOOST_TEST_MESSAGE( "test_terms_with_zero_descendants ... OK" );

}


BOOST_AUTO_TEST_CASE(test_term_descendants_BP)
{

  const std::string BP_term("GO:0051318");
  const std::set<std::string> descendants{"GO:0000080", "GO:0051330"};
  if (convertSet(goGraph().getGoGraphImpl().getDescendantTerms(BP_term)) != descendants) BOOST_FAIL( "BP term: " + BP_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_descendants_MF)
{

  const std::string MF_term("GO:0042165");
  const std::set<std::string> descendants{"GO:0031626", "GO:0042166"};
  if (convertSet(goGraph().getGoGraphImpl().getDescendantTerms(MF_term)) != descendants) BOOST_FAIL( "MF term: " + MF_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_descendants_CC)
{
  const std::string CC_term("GO:0030057");
  const std::set<std::string> descendants{"GO:0090635", "GO:0090636", "GO:0090637"};
  if (convertSet(goGraph().getGoGraphImpl().getDescendantTerms(CC_term)) != descendants) BOOST_FAIL( "CC term: " + CC_term + " descendants not found" );
  BOOST_TEST_MESSAGE( "test_term_descendants_CC ... OK" );

}

////////////////////////////////////////////////
// Parents, immediate ancestors
////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_term_parents_bad_id)
{

  const std::string Bad_term("GO:0022");
  if (not goGraph().getGoGraphImpl().getParentTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found parents" );
  BOOST_TEST_MESSAGE( "test_term_parents_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_BP)
{

  const std::string BP_term("GO:0022403");
  const std::set<std::string> parents{"GO:0044848"};
  if (convertSet(goGraph().getGoGraphImpl().getParentTerms(BP_term)) != parents) BOOST_FAIL( "BP term: " + BP_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_MF)
{

  const std::string MF_term("GO:0048037");
  const std::set<std::string> parents{"GO:0005488"};
  if (convertSet(goGraph().getGoGraphImpl().getParentTerms(MF_term)) != parents) BOOST_FAIL( "MF term: " + MF_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_parents_CC)
{

  const std::string CC_term("GO:0005911");
  const std::set<std::string> parents{"GO:0030054"};
  if (convertSet(goGraph().getGoGraphImpl().getParentTerms(CC_term)) != parents) BOOST_FAIL( "CC term: " + CC_term + " parents not found" );
  BOOST_TEST_MESSAGE( "test_term_parents_CC ... OK" );

}

////////////////////////////////////////////////
// Children, immediate descendants
////////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_term_children_bad_id)
{

  const std::string Bad_term("GO:00098");
  if (not goGraph().getGoGraphImpl().getChildTerms(Bad_term).empty()) BOOST_FAIL( "Bad ontology term: " + Bad_term + " found children" );
  BOOST_TEST_MESSAGE( "test_term_children_bad_id ... OK" );

}
BOOST_AUTO_TEST_CASE(test_term_children_BP)
{

  const std::string BP_term("GO:0051319");
  const std::set<std::string> children{"GO:0051331", "GO:0000085"};
  if (convertSet(goGraph().getGoGraphImpl().getChildTerms(BP_term)) != children) BOOST_FAIL( "BP term: " + BP_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_children_MF)
{
  const std::string MF_term("GO:0001071");
  const std::set<std::string> children{"GO:0003700", "GO:0001070"};
  if (convertSet(goGraph().getGoGraphImpl().getChildTerms(MF_term)) != children) BOOST_FAIL( "MF term: " + MF_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_term_children_CC)
{

  const std::string CC_term("GO:0031974");
  const std::set<std::string> children{"GO:1904724", "GO:1904813", "GO:0043233", "GO:0031970"};
  if (convertSet(goGraph().getGoGraphImpl().getChildTerms(CC_term)) != children) BOOST_FAIL( "CC term: " + CC_term + " children not found" );
  BOOST_TEST_MESSAGE( "test_term_children_CC ... OK" );

}

/////////////////////////////////////////////////////////////////////
// Retrieve terms sets and ontology term sets
/////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_all_terms_accessor)
{

  const size_t term_count = 42718;
  if (goGraph().getGoGraphImpl().getAllTerms().size() != term_count) BOOST_FAIL( "All terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_BP)
{

  const size_t term_count = 28651;
  if (goGraph().getGoGraphImpl().getAllTermsBP().size() != term_count) BOOST_FAIL( "All BP terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_MF)
{

  const size_t term_count = 10160;
  if (goGraph().getGoGraphImpl().getAllTermsMF().size() != term_count) BOOST_FAIL( "All MF terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_all_terms_accessor_CC)
{

  const size_t term_count = 3907;
  if (goGraph().getGoGraphImpl().getAllTermsCC().size() != term_count) BOOST_FAIL( "All CC terms count not equal : " + std::to_string(term_count));
  BOOST_TEST_MESSAGE( "test_all_terms_accessor_CC ... OK" );

}

///////////////////////////////////////////////////////////////////////////
// Test filtering functions
///////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_filter_for_BP)
{

  kol::OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0033631", "GO:0033627"};
  if (convertSet(goGraph().getGoGraphImpl().filterSetForBP(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all BP");
  BOOST_TEST_MESSAGE( "test_filter_for_BP ... OK" );

}


BOOST_AUTO_TEST_CASE(test_filter_for_MF)
{

  kol::OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0004099", "GO:0019213"};
  if (convertSet(goGraph().getGoGraphImpl().filterSetForMF(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all MF");
  BOOST_TEST_MESSAGE( "test_filter_for_MF ... OK" );

}


BOOST_AUTO_TEST_CASE(test_filter_for_CC)
{

  kol::OntologySetType<std::string> unfiltered_list{"GO:0033631", "GO:0033627", "GO:0004099", "GO:0019213", "GO:0044225", "GO:0043025"};
  std::set<std::string> filtered_list{"GO:0044225", "GO:0043025"};
  if (convertSet(goGraph().getGoGraphImpl().filterSetForCC(unfiltered_list)) != filtered_list) BOOST_FAIL( "Filtered terms count not all CC");
  BOOST_TEST_MESSAGE( "test_filter_for_CC ... OK" );

}

///////////////////////////////////////////////////////////
// Test accessors for the root terms of the 3 ontologies
//////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_root_term_access_BP)
{

  if (goGraph().getGoGraphImpl().getTermRoot("GO:0033631") != kol::GO::getRootTermBP()) BOOST_FAIL( "Failed to get root term for BP");
  BOOST_TEST_MESSAGE( "test_root_term_access_BP ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_term_access_MF)
{

  if (goGraph().getGoGraphImpl().getTermRoot("GO:0019213") != kol::GO::getRootTermMF()) BOOST_FAIL( "Failed to get root term for MF");
  BOOST_TEST_MESSAGE( "test_root_term_access_MF ... OK" );

}

BOOST_AUTO_TEST_CASE(test_root_term_access_CC)
{

  if (goGraph().getGoGraphImpl().getTermRoot("GO:0044225") != kol::GO::getRootTermCC()) BOOST_FAIL( "Failed to get root term for CC");
  BOOST_TEST_MESSAGE( "test_root_term_access_CC ... OK" );

}

BOOST_AUTO_TEST_SUITE_END()
