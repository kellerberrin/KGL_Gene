/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_JIANG_CONRATH_SIMILARITY
#define KGL_JIANG_CONRATH_SIMILARITY

#include "kol_TermSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include <kol_GoGraph.h>

namespace kellerberrin::ontology {


/*! \class JiangConrathSimilarity
	\brief A class to calculate Jiang Conrath similarity between 2 terms

	This class calculates Jiang Conrath similarity.
	
	Jiang, J. J., & Conrath, D. W. (1997). Semantic similarity based on corpus 
	  statistics and lexical taxonomy. In Proc. of 10th International Conference
	  on Research on Computational Linguistics, Taiwan.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	distance = IC(termA) + IC(termB) - 2*IC(MICA)
	maxDistance = 2*IC(single annotaiotn)
	similarity = 1 - distance/maxDistance
	(see Lord et al.)

*/
class JiangConrathSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  JiangConrathSimilarity( const std::shared_ptr<const GoGraph> &graph_ptr,
                          const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~JiangConrathSimilarity() override = default;
  //! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
  /*!
    This method returns the JiangConrath similarity or the information content of the most informative common ancestor.
  */

  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

  //! A method for calculating term-to-term similarity for GO terms using Normalized JiangConrath similarity
  /*!
    This method returns the JiangConrath similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //JiangConrath's method is already normalized
    return calculateTermSimilarity(goTermA, goTermB);

  }


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationContentMap> ic_map_ptr_;

};


} // namespace


#endif
