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


  [[nodiscard]] const SimilarityInterface &termSimilarityAnalysis() {

    if (not term_similarity_ptr_) {

      term_similarity_ptr_ = getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

protected:

  [[nodiscard]] std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] std::shared_ptr<const GoGraphImpl> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }



private:

  inline static std::unique_ptr<const SimilarityInterface> term_similarity_ptr_;
  inline const static TestValues test_values_;

  [[nodiscard]] virtual std::unique_ptr<const SimilarityInterface> getSimilarityAnalysis() = 0;


};


template<class SimAnalysis, class TestValues>
class TestTermSimilarity : public TestSimilarity<TestValues> {

public:

  TestTermSimilarity() = default;

  ~TestTermSimilarity() override = default;

private:

  [[nodiscard]] std::unique_ptr<const SimilarityInterface> getSimilarityAnalysis() override {

    std::shared_ptr<const GoGraphImpl> graph_ptr = TestSimilarity<TestValues>::getGoGraph();
    std::shared_ptr<const TermAnnotation> annotation_ptr = TestSimilarity<TestValues>::getAnnotation();
    std::shared_ptr<const InformationContentDAG> info_map_ptr(std::make_shared<const InformationContentDAG>(graph_ptr, annotation_ptr));
    return std::make_unique<const SimAnalysis>(info_map_ptr);

  }


};

template<class TestValues>
class TestPekarStaabSimilarity : public TestSimilarity<TestValues> {

public:

  TestPekarStaabSimilarity() = default;

  ~TestPekarStaabSimilarity() override = default;

private:

  [[nodiscard]] std::unique_ptr<const SimilarityInterface> getSimilarityAnalysis() override {

    std::shared_ptr<const GoGraphImpl> graph_ptr = TestSimilarity<TestValues>::getGoGraph();
    std::shared_ptr<const InformationDepthMap> info_map_ptr(std::make_shared<const InformationDepthMap>(*graph_ptr));
    return std::make_unique<const SimilarityPekarStaab>(graph_ptr, info_map_ptr);

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


  [[nodiscard]] const TermAnnotation &annotation() {

    if (not annotation_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(annotation_ptr_);
    return *annotation_ptr_;

  }

  [[nodiscard]] const SetSimilarityInterface &setSimilarity() {

    if (not set_similarity_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(set_similarity_ptr_);
    return *set_similarity_ptr_;

  }

  [[nodiscard]] const GoGraphImpl &goGraph() {

    if (not graph_ptr_) {

      getSetSimilarityAnalysis();

    }
    BOOST_REQUIRE(graph_ptr_);
    return *graph_ptr_;

  }


  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};


protected:

  inline static std::shared_ptr<const TermAnnotation> annotation_ptr_;
  inline static std::shared_ptr<const GoGraphImpl> graph_ptr_;
  inline static std::unique_ptr<const SetSimilarityInterface> set_similarity_ptr_;
  inline const static TestValues test_values_;

  virtual void getSetSimilarityAnalysis() = 0;

  [[nodiscard]] static std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] static std::shared_ptr<const TermAnnotation> getSymbolicAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile( default_evidence,
                                                     UnitTestDefinitions::gafFileName(),
                                                     AnnotationGeneName::SYMBOLIC_GENE_ID);

  }


  [[nodiscard]] static std::shared_ptr<const GoGraphImpl> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
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
    std::shared_ptr<const InformationContentDAG> ic_map_ptr(std::make_shared<const InformationContentDAG>(TSS::graph_ptr_, TSS::annotation_ptr_));
    std::shared_ptr<const SimAnalysis> similarity_ptr(std::make_shared<const SimAnalysis>(ic_map_ptr));
    TSS::set_similarity_ptr_ = std::make_unique<const SetSimAnalysis>(similarity_ptr);

  }


};


template<class SimAnalysis, class SetSimAnalysis, class TestValues>
class TestSetTermSymbolicSimilarity : public TestSetSimilarity<TestValues> {

public:

  TestSetTermSymbolicSimilarity() = default;

  ~TestSetTermSymbolicSimilarity() override = default;


private:

  using TSS = TestSetSimilarity<TestValues>;

  void getSetSimilarityAnalysis() override {

    TSS::graph_ptr_ = TSS::getGoGraph();
    TSS::annotation_ptr_ = TSS::getSymbolicAnnotation();
    std::shared_ptr<const InformationContentDAG> ic_map_ptr(std::make_shared<const InformationContentDAG>(TSS::graph_ptr_, TSS::annotation_ptr_));
    std::shared_ptr<const SimAnalysis> similarity_ptr(std::make_shared<const SimAnalysis>(ic_map_ptr));
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
    std::shared_ptr<const InformationContentDAG> info_map_ptr(std::make_shared<const InformationContentDAG>(TSS::graph_ptr_, TSS::annotation_ptr_));
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


  [[nodiscard]] const InformationInterface &sharedAnalysis() {

    if (not term_similarity_ptr_) {

      getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_similarity_ptr_);
    return *term_similarity_ptr_;

  }

  [[nodiscard]] const InformationContent &termInformation() {

    if (not term_information_ptr_) {

      getSimilarityAnalysis();

    }
    BOOST_REQUIRE(term_information_ptr_);
    return *term_information_ptr_;

  }


  [[nodiscard]] static const TestValues &getValue() { return test_values_; }

  const static constexpr double TEST_ACCURACY_PERCENT{0.0001};

protected:

  [[nodiscard]] std::shared_ptr<const TermAnnotation> getAnnotation() {

    PolicyEvidence default_evidence;
    return ParserAnnotationGaf::parseAnnotationFile(default_evidence, UnitTestDefinitions::gafFileName());

  }

  [[nodiscard]] std::shared_ptr<const GoGraphImpl> getGoGraph() {

    auto go_parser_ptr = ParserGoFactory::createGoParser(ParserGoType::PARSER_GO_OBO);
    BOOST_REQUIRE(go_parser_ptr);
    return go_parser_ptr->parseGoFile(UnitTestDefinitions::oboFileName());

  }

  inline static std::shared_ptr<const InformationContent> term_information_ptr_;
  inline static std::unique_ptr<const InformationInterface> term_similarity_ptr_;

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

    std::shared_ptr<const GoGraphImpl> graph_ptr = SS::getGoGraph();
    std::shared_ptr<const TermAnnotation> annotation_ptr = SS::getAnnotation();
    SS::term_information_ptr_ = std::make_shared<const InformationContent>(graph_ptr, annotation_ptr);
    SS::term_similarity_ptr_ = std::make_unique<const SharedAnalysis>(graph_ptr, SS::term_information_ptr_);

  }

};


} // namespace


#endif //GGTK_TEST_H
