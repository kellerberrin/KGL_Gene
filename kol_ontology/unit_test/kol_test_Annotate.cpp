//
// Created by kellerberrin on 3/4/21.
//

#include <kol_library.h>
#include "kol_test.h"
#include <boost/test/unit_test.hpp>

namespace kellerberrin::ontology {

// This object is re-created for each test case; store the graph in a static pointer so that it is only created once.
class TestDepthMap {

public:

  TestDepthMap() = default;

  ~TestDepthMap() = default;

  [[nodiscard]] static const TermAnnotation &annotation() {

    if (not static_annotation_) {

      static_annotation_ = getAnnotation();

    }
    BOOST_REQUIRE(static_annotation_);
    return *static_annotation_;

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

  [[nodiscard]] static std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }


  inline static std::shared_ptr<const TermAnnotation> static_annotation_;

};


} // namespace

namespace kol = kellerberrin::ontology;

BOOST_FIXTURE_TEST_SUITE(TestAnnotateSuite, kol::TestDepthMap)


///////////////////////////////////////////////
// Gene and GO Term count accessors
///////////////////////////////////////////////



BOOST_AUTO_TEST_CASE(test_gene_count_accessor)
{

  const size_t gene_count{19187};
  const size_t anno_gene_count = annotation().getNumGenes();
  if (anno_gene_count != gene_count) BOOST_FAIL("Annotation gene count is incorrect: " + std::to_string(anno_gene_count));
  BOOST_TEST_MESSAGE( "test_gene_count_accessor... OK" );

}

BOOST_AUTO_TEST_CASE(test_go_count_accessor)
{

  const size_t go_count{16496};
  const size_t anno_go_terms = annotation().getNumGoTerms();
  if (anno_go_terms != go_count) BOOST_FAIL("Annotation GO count is incorrect: " + std::to_string(anno_go_terms));
  BOOST_TEST_MESSAGE( "test_go_count_accessor ... OK" );

}

///////////////////////////////////////////////////
// Access all Genes and GO Term as lists
//////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_get_all_go_terms)
{

        const size_t go_count{16496};
        // Check unqiueness by assigning to a set.
        auto go_set = annotation().getAllGoTerms();
        if (go_set.size() != go_count) BOOST_FAIL("Annotation Go list size is incorrect: " + std::to_string(go_set.size()));
        BOOST_TEST_MESSAGE( "test_get_all_go_terms... OK" );

}

BOOST_AUTO_TEST_CASE(test_get_all_genes)
{

  const size_t gene_count{19187};
  // Check unqiueness by assigning to a set.
  auto gene_set = annotation().getAllGenes();
  if (gene_set.size() != gene_count) BOOST_FAIL("Annotation Gene list size is incorrect: " + std::to_string(gene_set.size()));
  BOOST_TEST_MESSAGE( "test_get_all_genes ... OK" );

}

//////////////////////////////////////////////////////////
// Test existence of genes and go terms in the database
//////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(test_has_gene_bad_id)
{

  if (annotation().hasGene("ABC12345")) BOOST_FAIL("Annotation bad gene id found" );
  BOOST_TEST_MESSAGE( "test_has_gene_bad_id... OK" );

}

BOOST_AUTO_TEST_CASE(test_has_gene)
{

  if (not annotation().hasGene("A0A024R161")) BOOST_FAIL("Annotation valid gene id found" );
  BOOST_TEST_MESSAGE( "test_has_gene ... OK" );

}

////////////////////////////////////////////////
// Per Gene/ Per Term annotation counts
////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_num_gene_annotations_bad_id)
{

  if (annotation().getNumAnnotationsForGene("ABC12345") != 0) BOOST_FAIL("Found gene annotations for bad gene id" );
  BOOST_TEST_MESSAGE( "test_num_gene_annotations_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_num_gene_annotations)
{

  const size_t gene_annotations{3};
  if (annotation().getNumAnnotationsForGene("A0A024R161") != gene_annotations) BOOST_FAIL("Found incorrect number of gene annotations." );
  BOOST_TEST_MESSAGE( "test_num_gene_annotations ... OK" );

}

BOOST_AUTO_TEST_CASE(test_num_go_annotations_bad_id)
{

  if (annotation().getNumAnnotationsForGoTerm("GO:00") != 0) BOOST_FAIL("Found go annotations for bad gene id" );
  BOOST_TEST_MESSAGE( "test_num_go_annotations_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_num_go_annotations)
{

//  const size_t gaf_record_annotations{216};
  const size_t go_annotations{196};
  auto gene_map =  annotation().getGOGeneMap("GO:0004871");
  size_t gaf_count{0};
  for (auto const& [gene_id, gaf_vector] : gene_map) {

    gaf_count += gaf_vector.size();

  }
  size_t found_annotations = annotation().getNumAnnotationsForGoTerm("GO:0004871");
  if (found_annotations != go_annotations) BOOST_FAIL("Found incorrect number of go annotations: "
                                                      + std::to_string(found_annotations)
                                                      + "  GAF Count: " + std::to_string(gaf_count));
  BOOST_TEST_MESSAGE( "test_num_go_annotations ... OK" );

}

//////////////////////////////////////////////////////
// List of go terms for genes and vice versa
//////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(test_go_terms_for_gene_bad_id)
{

  if (not annotation().getGoTermsForGene("ABC12345").empty()) BOOST_FAIL("Found go annotation for bad gene id" );
  BOOST_TEST_MESSAGE( "test_go_terms_for_gene_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_go_terms_for_gene)
{

  const std::set<std::string> go_annotations{"GO:0004871", "GO:0005834", "GO:0007186"};
  if (convertSet(kol::SetUtilities::convertVector(annotation().getGoTermsForGene("A0A024R161"))) != go_annotations) BOOST_FAIL("Found invalid go terms for gene" );
  BOOST_TEST_MESSAGE( "test_go_terms_for_gene ... OK" );

}

BOOST_AUTO_TEST_CASE(test_genes_for_go_term_bad_id)
{

  if (not annotation().getGenesForGoTerm("ABC12345").empty()) BOOST_FAIL("Found genes for bad go term" );
  BOOST_TEST_MESSAGE( "test_genes_for_go_term_bad_id ... OK" );

}

BOOST_AUTO_TEST_CASE(test_genes_for_go_term)
{

  const std::set<std::string> gene_set{"A6NNW6", "P06733", "P09104", "P13929"};
  if (convertSet(kol::SetUtilities::convertVector(annotation().getGenesForGoTerm("GO:0000015"))) != gene_set) BOOST_FAIL("Invalid gene list for go term" );
  BOOST_TEST_MESSAGE( "test_genes_for_go_term ... OK" );

}




BOOST_AUTO_TEST_CASE(test_evidence_codes_and_go_list_w_gene_query)
{

  const kol::OntologySetType<std::string> go_terms = {"GO:0005813", "GO:0007224", "GO:0007286",
                                             "GO:0030992", "GO:0036064", "GO:0042073",
                                             "GO:0042384", "GO:0072372", "GO:0072372",
                                             "GO:0072372", "GO:0072372", "GO:0072372",
                                             "GO:0097542", "GO:0097542", "GO:0097542", "GO:0097542"};


  auto gene_go_terms = annotation().getGoTermsForGene("A0AVF1");
  auto set_go_terms = kol::SetUtilities::convertVector(gene_go_terms);
  if (set_go_terms != go_terms) {

    BOOST_FAIL("Found invalid go terms for gene" );

  }
  BOOST_TEST_MESSAGE( "test_evidence_codes_and_go_list_w_gene_query ... OK" );

}


BOOST_AUTO_TEST_SUITE_END()
