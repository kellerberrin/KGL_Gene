//
// Created by kellerberrin on 7/4/21.
//

#ifndef KGL_KOL_TERMSIMILARITYCACHE_H
#define KGL_KOL_TERMSIMILARITYCACHE_H


#include <vector>
#include <string>

#include "kol_GoEnums.h"
#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermSimilarityInterface.h"

namespace kellerberrin::ontology {


//! A multi-threaded class write a term similarity matrix to a memory cache.
/*! \class TermSimilarityCache
	This class creates a memory cache similarity matrix by calculating all similarity values
	  between pairs of terms.
	  Warning! this object can use several gigabytes of memory.
	  Matrix creation time will depend on the number of execution threads committed to matrix creation.
	  Defaults to (HW threads available - 1).
*/
class TermSimilarityCache : public TermSimilarityInterface {


public:
  //! Parameterized constructor
  /*!
    A simple parameterized constructor.
    This class takes an instance of a similarity interface.
  */
  TermSimilarityCache(const std::shared_ptr<const GoGraph> &go_graph_ptr,
                      const std::shared_ptr<const AnnotationData> &annotation_ptr,
                      const std::shared_ptr<const TermSimilarityInterface> &term_sim_ptr,
                      GO::Ontology ontology = GO::Ontology::BIOLOGICAL_PROCESS) {

    termSimilarityCache(go_graph_ptr, annotation_ptr, term_sim_ptr, ontology);

  }

  TermSimilarityCache(const TermSimilarityCache &) = delete; // Too expensive.
  ~TermSimilarityCache() override = default;

  TermSimilarityCache &operator=(const TermSimilarityCache &) = delete; // Too expensive.

  //! A method for calculating term-to-term similarity for GO terms using a precomputed similarity matrix.
  /*!
    This method returns the term similarity as defined by the matrix.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &row_term, const std::string &column_term) const override;


  // termCount() == 0 is an error condition.
  [[nodiscard]] size_t termCount() const { return cache_matrix_.size(); }

private:

  OntologyMapType<std::string, std::size_t> term_to_index_;
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
  bool termSimilarityCache(const std::shared_ptr<const GoGraph> &go_graph_ptr,
                           const std::shared_ptr<const AnnotationData> &annotation_ptr,
                           const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                           GO::Ontology ontology);

  // Multithreaded calculation.
  static std::unique_ptr<std::vector<double>> calcColumn(size_t row,
                                                         const std::string &term_A,
                                                         const std::shared_ptr<const TermSimilarityInterface> &term_similarity_ptr,
                                                         const std::shared_ptr<const std::vector<std::string>> &ontology_terms_ptr);
};


} // namespace


#endif //KGL_TERMSIMILARITYCACHE_H
