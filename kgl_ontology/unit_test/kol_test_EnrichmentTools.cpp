//
// Created by kellerberrin on 3/4/21.
//

#include <kol_library.h>
#include "kol_test.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestEnrichment {

public:

  TestEnrichment() = default;

  ~TestEnrichment() = default;

  [[nodiscard]] static const AnnotationData &annotation() {

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

  inline static std::unique_ptr<const GoGraph> static_graph_;
  inline static std::unique_ptr<const AnnotationData> static_annotation_;

  [[nodiscard]] static std::unique_ptr<AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::unique_ptr<GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
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
  auto child_terms = goGraph().getDescendantTerms("GO:0002507");
  child_terms.insert("GO:0002507");

  for (auto const& gene : term_genes) {

    auto go_BP = annotation().getGoTermsForGeneBP(gene, goGraph());
    if (SetUtilities::set_intersection(go_BP, child_terms).empty()) BOOST_FAIL( "Empty intersection of BP go terms and child go terms" );

  }
  BOOST_TEST_MESSAGE( "test_enrichment_tools_child_genes ... OK" );

}

BOOST_AUTO_TEST_CASE(test_enrichment_tools_enrichment_significance)
{

  auto term_genes = kol::EnrichmentTools::getDescendantGenes(goGraph(), annotation(), "GO:0002507");
  auto p_val = kol::EnrichmentTools::enrichmentSignificance(goGraph(), annotation(), term_genes, "GO:0002507");
  BOOST_CHECK_CLOSE( p_val, 1.9224043e-43, 0.0001);
  BOOST_TEST_MESSAGE( "test_enrichment_tools_enrichment_significance ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
