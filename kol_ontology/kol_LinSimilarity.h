/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_LIN_SIMILARITY
#define KGL_LIN_SIMILARITY

#include "kol_TermSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {

/*! \class linSimPtr
	\brief A class to calculate Lin similarity between 2 terms

	This class calculates Lin similarity.
	
	Lin, D. (1998) An information theoretic definition of similarity. In: Proc. of the
	  15th International Conference on Machine Learning. San Francisco, CA:
	  Morgan Kaufman. pp 296-304

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	2 * IC(MICA) / ( IC(termA) + IC(termB) )

*/
class LinSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the lin Similarity Calculator.
  */
  LinSimilarity( const std::shared_ptr<const GoGraph> &graph_ptr,
                 const std::shared_ptr<const TermInformationInterface> &ic_map_ptr_)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr_) {}

  ~LinSimilarity() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Lin similarity
  /*!
    This method returns the Lin similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

  //! A method for calculating term-to-term similarity for GO terms using Normalized Lin similarity
  /*!
    This method returns the Lin similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //Lin's method is already normalized
    return calculateTermSimilarity(goTermA, goTermB);

  }

private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationInterface> ic_map_ptr_;

};

} // namespace

#endif
