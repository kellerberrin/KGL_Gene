//
// Created by kellerberrin on 2/4/21.
//

#ifndef GGTK_TEST_H
#define GGTK_TEST_H

#include "ggtk.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the necessary file and directory definitions.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
class UnitTestDefinitions {

public:

  // Just static definitions.
  UnitTestDefinitions() = delete;

  // GoGraph file names.
  [[nodiscard]] static std::string oboFileName() { return std::string(GRAPH_DIRECTORY_) + std::string(GO_OBO_); }
  [[nodiscard]] static std::string xmlFileName() { return std::string(GRAPH_DIRECTORY_) + std::string(GO_XML_); }
  // Association file names.
  [[nodiscard]] static std::string gafFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(HUMAN_GAF_ANNOTATION_); }
  [[nodiscard]] static std::string geneFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(GENE_ASSOCIATION_); }
  [[nodiscard]] static std::string entrezFileName() { return std::string(ANNOTATION_DIRECTORY_) + std::string(HUMAN_ENTREZ_ANNOTATION_); }

private:

  // Modify these directory constants to change the relative (or absolute) location of the data files.
  static const constexpr char* GRAPH_DIRECTORY_ = "Additional/ggtk/example_graphs/";
  static const constexpr char* ANNOTATION_DIRECTORY_ = "Additional/ggtk/example_annotations/";


  static const constexpr char* GO_OBO_ = "go-basic.obo";
  static const constexpr char* GO_XML_ = "go_daily-termdb.obo-xml";

  static const constexpr char* HUMAN_GAF_ANNOTATION_ = "goa_human.gaf";
  static const constexpr char* GENE_ASSOCIATION_ = "gene_association.mgi";
  static const constexpr char* HUMAN_ENTREZ_ANNOTATION_ = "human_gene2go";

};

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the similarity analysis object in a static pointer so that it is only created once.
/////////////////////////////////////////////////////////////////////////////////////////
template<class TestValues>
class TestSimilarity
{

public:

  TestSimilarity() = default;
  virtual ~TestSimilarity() = default;


  [[nodiscard]] const TermSimilarityInterface& termSimilarityAnalysis() {

    if (not term_similarity_ptr_) {

      term_similarity_ptr_ = getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }
  [[nodiscard]] static const TestValues& getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

protected:

  [[nodiscard]] std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

private:

  inline static std::unique_ptr<const TermSimilarityInterface> term_similarity_ptr_;
  const static TestValues test_values_;

  [[nodiscard]] virtual std::unique_ptr<const TermSimilarityInterface> getSimilarityAnalysis() = 0;


};



template<class SimAnalysis, class TestValues>
class TestTermSimilarity : public TestSimilarity<TestValues> {

public:

  TestTermSimilarity() = default;
  ~TestTermSimilarity() override = default;

private:

  [[nodiscard]] std::unique_ptr<const TermSimilarityInterface> getSimilarityAnalysis() override {

    std::shared_ptr<const GoGraph> graph_ptr = TestSimilarity<TestValues>::getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = TestSimilarity<TestValues>::getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(graph_ptr, annotation_ptr));
    return std::make_unique<const SimAnalysis>(graph_ptr, info_map_ptr);

  }

};

template<class TestValues>
class TestPekarStaabSimilarity : public TestSimilarity<TestValues> {

public:

  TestPekarStaabSimilarity() = default;
  ~TestPekarStaabSimilarity() override = default;

private:

  [[nodiscard]] std::unique_ptr<const TermSimilarityInterface> getSimilarityAnalysis() override {

    std::shared_ptr<const GoGraph> graph_ptr = TestSimilarity<TestValues>::getGoGraph();
    std::shared_ptr<const TermDepthMap> info_map_ptr(std::make_shared<const TermDepthMap>(*graph_ptr));
    return std::make_unique<const PekarStaabSimilarity>(graph_ptr, info_map_ptr);

  }


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test Modular Similarity Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template<class SimAnalysis, class  ModSimAnalysis>
class TestModTermSimilarity
{

public:

  TestModTermSimilarity() = default;
  ~TestModTermSimilarity() = default;


  [[nodiscard]] static const TermSimilarityInterface& termSimilarity() {

    if (not term_similarity_ptr_) {

      getModSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] static const TermSimilarityInterface& termModSimilarity() {

    if (not term_modular_similarity_ptr_) {

      getModSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_modular_similarity_ptr_);
    return *term_modular_similarity_ptr_;

  }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};


private:

  inline static std::unique_ptr<const TermSimilarityInterface> term_similarity_ptr_;
  inline static std::unique_ptr<const TermSimilarityInterface> term_modular_similarity_ptr_;

  static void getModSimilarityAnalysis() {

    std::shared_ptr<const GoGraph> graph_ptr = getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(graph_ptr, annotation_ptr));
    term_similarity_ptr_ = std::make_unique<const SimAnalysis>(graph_ptr, info_map_ptr);
    std::shared_ptr<const MICASharedInformation> shared_information_ptr(std::make_shared<const MICASharedInformation>(graph_ptr, info_map_ptr));
    term_modular_similarity_ptr_ = std::make_unique<const ModSimAnalysis>(shared_information_ptr);

  }

  [[nodiscard]] static std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test Set Similarity Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class SimAnalysis, class  SetSimAnalysis, class TestValues>
class TestSetTermSimilarity
{

public:

  TestSetTermSimilarity() = default;
  ~TestSetTermSimilarity() = default;


  [[nodiscard]] static const AnnotationData& annotation() {

    if (not annotation_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(annotation_ptr_);
    return *annotation_ptr_;

  }

  [[nodiscard]] static const TermSetSimilarityInterface& setSimilarity() {

    if (not set_similarity_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(set_similarity_ptr_);
    return *set_similarity_ptr_;

  }

  [[nodiscard]] static const GoGraph& goGraph() {

    if (not graph_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(graph_ptr_);
    return *graph_ptr_;

  }


  [[nodiscard]] static const TestValues& getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};


private:

  inline static std::shared_ptr<const AnnotationData> annotation_ptr_;
  inline static std::shared_ptr<const GoGraph> graph_ptr_;
  inline static std::unique_ptr<const SetSimAnalysis> set_similarity_ptr_;
  const static TestValues test_values_;

  static void getSetSimilarityAnalysis() {

    graph_ptr_ = getGoGraph();
    annotation_ptr_ = getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(graph_ptr_, annotation_ptr_));
    std::shared_ptr<const SimAnalysis> similarity_ptr(std::make_shared<const SimAnalysis>(graph_ptr_, info_map_ptr));
    set_similarity_ptr_ = std::make_unique<const SetSimAnalysis>(similarity_ptr);

  }

  [[nodiscard]] static std::unique_ptr<const AnnotationData> getAnnotation() {

    auto anno_parser_ptr = AnnotationParserFactory::createAnnotationParser(AnnotationParserType::GAF_ANNO_PARSER,
                                                                           DisallowedSetEvidencePolicy());
    BOOST_REQUIRE(anno_parser_ptr);
    return anno_parser_ptr->parseAnnotationFile(UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::unique_ptr<const GoGraph> getGoGraph() {

    auto go_parser_ptr = GoParserFactory::createGoParser(GoParserType::OBO_GO_STANDARD);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

};



#endif //GGTK_TEST_H
