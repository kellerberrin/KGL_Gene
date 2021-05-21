//
// Created by kellerberrin on 2/4/21.
//

#ifndef KGL_KOL_TEST_H
#define KGL_KOL_TEST_H

#include "kol_library.h"
#include "kol_test_data.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>


namespace kellerberrin::ontology {

///////////////////////////////////////////////////////////////////////////////////////
// This object is re-created for each test case.
// Store the similarity analysis object in a static pointer so that it is only created once.
/////////////////////////////////////////////////////////////////////////////////////////
template<class TestValues>
class TestSimilarity {

public:

  TestSimilarity() = default;

  virtual ~TestSimilarity() = default;


  [[nodiscard]] const TermSimilarityInterface &termSimilarityAnalysis() {

    if (not term_similarity_ptr_) {

      term_similarity_ptr_ = getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

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
  inline const static TestValues test_values_;

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
    std::shared_ptr<const MICASharedInformation> shared_information_ptr(std::make_shared<const MICASharedInformation>(graph_ptr, info_map_ptr));
    return std::make_unique<const SimAnalysis>(shared_information_ptr);

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




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test Set Similarity Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TestValues>
class TestSetSimilarity {

public:

  TestSetSimilarity() = default;

  virtual ~TestSetSimilarity() = default;


  [[nodiscard]] const AnnotationData &annotation() {

    if (not annotation_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(annotation_ptr_);
    return *annotation_ptr_;

  }

  [[nodiscard]] const TermSetSimilarityInterface &setSimilarity() {

    if (not set_similarity_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(set_similarity_ptr_);
    return *set_similarity_ptr_;

  }

  [[nodiscard]] const GoGraph &goGraph() {

    if (not graph_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(graph_ptr_);
    return *graph_ptr_;

  }


  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};


protected:

  inline static std::shared_ptr<const AnnotationData> annotation_ptr_;
  inline static std::shared_ptr<const GoGraph> graph_ptr_;
  inline static std::unique_ptr<const TermSetSimilarityInterface> set_similarity_ptr_;
  inline const static TestValues test_values_;

  virtual void getSetSimilarityAnalysis() = 0;

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


template<class SimAnalysis, class SetSimAnalysis, class TestValues>
class TestSetTermSimilarity : public TestSetSimilarity<TestValues> {

public:

  TestSetTermSimilarity() = default;

  ~TestSetTermSimilarity() override = default;


private:

  using TSS = TestSetSimilarity<TestValues>;

  void getSetSimilarityAnalysis() override {

    TSS::graph_ptr_ = TSS::getGoGraph();
    TSS::annotation_ptr_ = TSS::getAnnotation();
    std::shared_ptr<const TermInformationContentMap> ic_map_ptr(std::make_shared<const TermInformationContentMap>(TSS::graph_ptr_, TSS::annotation_ptr_));
    std::shared_ptr<const MICASharedInformation> info_map_ptr(std::make_shared<const MICASharedInformation>( TSS::graph_ptr_, ic_map_ptr));
    std::shared_ptr<const SimAnalysis> similarity_ptr(std::make_shared<const SimAnalysis>(info_map_ptr));
    TSS::set_similarity_ptr_ = std::make_unique<const SetSimAnalysis>(similarity_ptr);

  }


};


template<class TestValues>
class TestJaccardSetSimilarity : public TestSetSimilarity<TestValues> {

public:

  TestJaccardSetSimilarity() = default;

  ~TestJaccardSetSimilarity() override = default;


private:

  using TSS = TestSetSimilarity<TestValues>;

  void getSetSimilarityAnalysis() override {

    TSS::graph_ptr_ = TSS::getGoGraph();
    TSS::annotation_ptr_ = TSS::getAnnotation();
    TSS::set_similarity_ptr_ = std::make_unique<const JaccardSetSimilarity>();

  }


};


template<class SetSimilarity, class TestValues>
class TestGraphSetSimilarity : public TestSetSimilarity<TestValues> {

public:

  TestGraphSetSimilarity() = default;

  ~TestGraphSetSimilarity() override = default;


private:

  using TSS = TestSetSimilarity<TestValues>;

  void getSetSimilarityAnalysis() override {

    TSS::graph_ptr_ = TSS::getGoGraph();
    TSS::annotation_ptr_ = TSS::getAnnotation();
    TSS::set_similarity_ptr_ = std::make_unique<const SetSimilarity>(TSS::graph_ptr_);

  }


};


template<class SetSimilarity, class TestValues>
class TestSetTheorySimilarity : public TestSetSimilarity<TestValues> {

public:

  TestSetTheorySimilarity() = default;

  ~TestSetTheorySimilarity() override = default;


private:

  using TSS = TestSetSimilarity<TestValues>;

  void getSetSimilarityAnalysis() override {

    TSS::graph_ptr_ = TSS::getGoGraph();
    TSS::annotation_ptr_ = TSS::getAnnotation();
    std::shared_ptr<const TermInformationContentMap> info_map_ptr(std::make_shared<const TermInformationContentMap>(TSS::graph_ptr_, TSS::annotation_ptr_));
    TSS::set_similarity_ptr_ = std::make_unique<const SetSimilarity>(TSS::graph_ptr_, info_map_ptr);

  }


};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Test Shared information Types
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<class TestValues>
class SharedSimilarity {

public:

  SharedSimilarity() = default;

  virtual ~SharedSimilarity() = default;


  [[nodiscard]] const SharedInformationInterface &sharedAnalysis() {

    if (not term_similarity_ptr_) {

      getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] const TermInformationContentMap &termInformation() {

    if (not term_information_ptr_) {

      getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_information_ptr_);
    return *term_information_ptr_;

  }


  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

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

  inline static std::shared_ptr<const TermInformationContentMap> term_information_ptr_;
  inline static std::unique_ptr<const SharedInformationInterface> term_similarity_ptr_;

private:

  inline const static TestValues test_values_;

  virtual void getSimilarityAnalysis() = 0;


};


template<class SharedAnalysis, class TestValues>
class TestSharedSimilarity : public SharedSimilarity<TestValues> {

public:

  TestSharedSimilarity() = default;

  ~TestSharedSimilarity() override = default;

private:

  using SS = SharedSimilarity<TestValues>;

  void getSimilarityAnalysis() override {

    std::shared_ptr<const GoGraph> graph_ptr = SS::getGoGraph();
    std::shared_ptr<const AnnotationData> annotation_ptr = SS::getAnnotation();
    SS::term_information_ptr_ = std::make_shared<const TermInformationContentMap>(graph_ptr, annotation_ptr);
    SS::term_similarity_ptr_ = std::make_unique<const SharedAnalysis>(graph_ptr, SS::term_information_ptr_);

  }

};


} // namespace


#endif //GGTK_TEST_H
