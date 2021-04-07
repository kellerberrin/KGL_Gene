//
// Created by kellerberrin on 7/4/21.
//

#ifndef TERMSIMILARITYCACHE_H
#define TERMSIMILARITYCACHE_H


#include <vector>
#include <string>
#include <iostream>

#include <GoEnums.h>
#include <GoGraph.h>
#include <AnnotationData.h>
#include <TermSimilarityInterface.h>

//! A multi-threaded class write a term similarity matrix to a memory cache.
/*! \class TermSimilarityCache
	This class creates a memory cache similarity matrix by calculating all similarity values
	  between pairs of terms. Warning this object can use several gigabytes of memory.
	  Matrix creation time will depend on the number of execution threads committed to matrix creation.
	  Defaults to (HW threads available - 1).
*/
class TermSimilarityCache : public TermSimilarityInterface {


public:
  //! Parameterized constructor
  /*!
    A simple parameterized constructor.
    This class take an instance of the GO Graph.
  */
  TermSimilarityCache(const std::shared_ptr<const TermSimilarityInterface>& term_sim_ptr,
                      GO::Ontology ontology = GO::Ontology::BIOLOGICAL_PROCESS) {

    cacheSimilarityMatrix(term_sim_ptr, ontology);

  }
  ~TermSimilarityCache() = default;


  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the term similarity as defined by the matrix.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

    if (hasTerm(goTermA) && hasTerm(goTermB)){

      auto const& [aterm, row] = *(_termToIndex.find(goTermA));
      auto const& [bterm, column] = *(_termToIndex.find(goTermB));

      return _matrix[row][column];

    } else {

      return 0.0;

    }

  }

  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

    return calculateTermSimilarity(goTermA,goTermB);

  }

  // termCount() == 0 is an error condition.
  [[nodiscard]] size_t termCount() const { return _matrix.size(); }

private:

  [[nodiscard]] bool hasTerm(const std::string &term) const {

    return _termToIndex.find(term) != _termToIndex.end();

  }

  //! A method to write a term similarity matrix
  /*!
    Calculates and writes the similarity matrix to file.
    Calculates the similarity between all pairs of terms.
    The complexity is O(N^2) * O(Term Pair Calculation Cost).
    O(Term Pair Calculation Cost) is usually near constant time,
     but some method will be extremely slow and require other methods
  */
  bool cacheSimilarityMatrix(const std::shared_ptr<const TermSimilarityInterface>& termSim,
                             GO::Ontology ontology) {

    // Return the empty cache.
    if (ontology == GO::Ontology::ONTO_ERROR) {

      return false;

    }

    std::vector<std::string> ontologyTerms = _annoData->getOntologyTerms(*_goGraph, ontology);

    size_t index{0};
    for (auto const& term : ontologyTerms) {

      auto const [iter, result] = _termToIndex.try_emplace(term, index);
      // Check that all terms are distinct
      if (not result) {
        // Return the empty cache.
        _termToIndex.clear();
        _matrix.clear();
        return false;

      }
      ++index;

    }

    // Create a cache matrix.
    size_t nTerms = ontologyTerms.size();
    _matrix.reserve(nTerms);

    for (size_t i = 0; i < nTerms; ++i){

      std::vector<double> column_vector(nTerms, 0.0);;
      _matrix.emplace_back(std::move(column_vector));

    }

    // --- Main Calculation ---
    for (size_t i = 0; i < nTerms; ++i){

      std::string termA = ontologyTerms.at(i);
      _matrix[i][i] = termSim->calculateNormalizedTermSimilarity(termA, termA);

      for (std::size_t j = i + 1; j < nTerms; ++j){

        std::string termB = ontologyTerms.at(j);
        double value = termSim->calculateNormalizedTermSimilarity(termA, termB);
        //cout << termA << "\t" << termB << "\t" << value << endl;
        _matrix[i][j] = value;
        _matrix[j][i] = value;

      }

    }

    return true;

  }

  OntologyMapType<std::string, std::size_t> _termToIndex;
  std::vector<std::vector<double> > _matrix;
  //! A pointer to the gene ontology graph
  std::shared_ptr<const GoGraph> _goGraph;
  //! A pointer to the annotation database
  std::shared_ptr<const AnnotationData> _annoData;

};




#endif //TERMSIMILARITYCACHE_H
