//
// Created by kellerberrin on 17/4/21.
//

#ifndef KOL_ASYMMETRICSIMILARITYCACHE_H
#define KOL_ASYMMETRICSIMILARITYCACHE_H


#include <vector>
#include <string>

#include "kol_GoEnums.h"
#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermSimilarityInterface.h"

namespace kellerberrin::ontology {


//! A multi-threaded class write an asymmetric term similarity matrix to a memory cache.
/*! \class AsymmetricSimilarityCache
	This class creates a memory cache similarity matrix by calculating all similarity values
	  between pairs of terms. Warning this object can use several gigabytes of memory.
	  Matrix creation time will depend on the number of execution threads committed to matrix creation.
	  Defaults to (HW threads available - 1).
*/
class AsymmetricSimilarityCache : public TermSimilarityInterface {


public:
  //! Parameterized constructor
  /*!
    A simple parameterized constructor.
    This class takes an instance of a similarity interface.
  */
  AsymmetricSimilarityCache(const std::vector<std::string> &row_terms,
                            const std::vector<std::string> &column_terms,
                            const std::shared_ptr<const TermSimilarityInterface> &term_sim_ptr,
                            GO::Ontology ontology = GO::Ontology::BIOLOGICAL_PROCESS) {

    termSimilarityCache(row_terms, column_terms, term_sim_ptr, ontology);

  }

  AsymmetricSimilarityCache(const AsymmetricSimilarityCache &) = delete; // Too expensive.
  ~AsymmetricSimilarityCache() override = default;

  AsymmetricSimilarityCache &operator=(const AsymmetricSimilarityCache &) = delete; // Too expensive.

  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the term similarity as defined by the matrix.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &row_term, const std::string &column_term) const override;


  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &row_term, const std::string &column_term) const override {

    return calculateTermSimilarity(row_term, column_term);

  }

  // rows() == 0 is an error condition.
  [[nodiscard]] size_t rows() const { return row_term_to_index_.size(); }
  [[nodiscard]] size_t columns() const { return column_term_to_index_.size(); }

private:

  OntologyMapType<std::string, std::size_t> row_term_to_index_;
  OntologyMapType<std::string, std::size_t> column_term_to_index_;
  std::vector<std::unique_ptr<std::vector<double>>> cache_matrix_;

  // Lookup function knows the matrix is symmetric
  [[nodiscard]] double lookUpMatrix(size_t row, size_t column) const;

  //! A method to write a term similarity cache
  /*!
    Calculates square matrix of similarity terms and caches this in memory.
    Calculates the similarity between all pairs of terms.
    The complexity is O(N^2) * O(Term Pair Calculation Cost).
    O(Term Pair Calculation Cost) is usually near constant time,
     but some methods will be extremely slow and infeasible.
  */
  bool termSimilarityCache(const std::vector<std::string> &row_terms,
                           const std::vector<std::string> &column_terms,
                           const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                           GO::Ontology ontology);

  // Multithreaded calculation.
  static std::unique_ptr<std::vector<double>> calcColumn(const std::string &row_term,
                                                         const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                         const std::shared_ptr<const std::vector<std::string>>& column_terms_ptr);
};


} // namespace


#endif //KOL_ASYMMETRICSIMILARITYCACHE_H
