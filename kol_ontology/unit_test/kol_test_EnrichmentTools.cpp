//
// Created by kellerberrin on 3/4/21.
//

#include <kol_library.h>
#include "kol_test.h"
#include "contrib/kol_GoGraphImpl.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestEnrichment {

public:

  TestEnrichment() = default;

  ~TestEnrichment() = default;

  [[nodiscard]] static const TermAnnotation &annotation() {

    if (not static_annotation_) {

      static_annotation_ = getAnnotation();

    }
    BOOST_REQUIRE(static_annotation_);
    return *static_annotation_;

  }

  [[nodiscard]] static const GoGraph &goGraph() {

    if (not static_graph_) {

      static_graph_ = getGoGraph();

    }
    BOOST_REQUIRE(static_graph_);
    return *static_graph_;

  }


private:

  inline static std::shared_ptr<const GoGraph> static_graph_;
  inline static std::shared_ptr<const TermAnnotation> static_annotation_;

  [[nodiscard]] static std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::shared_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }


};

} // namespace

namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;

BOOST_FIXTURE_TEST_SUITE(TestEnrichmentSuite, kol::TestEnrichment)


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_enrichment_tools_hypergeometric_test)
{

  kel::HypergeometricDistribution  hyper_dist(8,10,17);
  auto p_val = hyper_dist.upperSingleTailTest(7);
  BOOST_CHECK_CLOSE( p_val, 0.036404771, 0.0001);
  BOOST_TEST_MESSAGE( "test_enrichment_tools_hypergeometric_test ... OK" );

}

BOOST_AUTO_TEST_CASE(test_enrichment_tools_child_genes)
{

  auto term_genes = kol::EnrichmentTools::getDescendantGenes(goGraph(), annotation(), "GO:0002507");
  auto child_terms = goGraph().getGoGraphImpl().getDescendantTerms("GO:0002507");
  child_terms.insert("GO:0002507");

  for (auto const& gene : term_genes) {

    auto go_BP = annotation().getGoTermsForGeneBP(gene);
    if (kol::SetUtilities::setIntersection(go_BP, child_terms).empty()) BOOST_FAIL("Empty intersection of BP go terms and child go terms" );

  }
  BOOST_TEST_MESSAGE( "test_enrichment_tools_child_genes ... OK" );

}

BOOST_AUTO_TEST_CASE(test_enrichment_tools_enrichment_significance)
{

  auto term_genes = kol::EnrichmentTools::getDescendantGenes(goGraph(), annotation(), "GO:0002507");
  auto p_val = kol::EnrichmentTools::enrichmentSignificance(goGraph(), annotation(), term_genes, "GO:0002507");
  BOOST_CHECK_CLOSE( p_val, 0.0, 0.0001);
  BOOST_TEST_MESSAGE( "test_enrichment_tools_enrichment_significance ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
