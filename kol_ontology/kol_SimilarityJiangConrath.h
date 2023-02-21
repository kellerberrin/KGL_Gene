/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_MODULAR_JIANG_CONRATH
#define KOL_MODULAR_JIANG_CONRATH

#include "kol_SimilarityInterface.h"
#include "kol_InformationInterface.h"

#include <memory>

namespace kellerberrin::ontology {


/*! \class SimilarityJiangConrath
	\brief A class to calculate Jiang Conrath similarity between 2 terms

	This class calculates Jiang Conrath similarity.
	
	Jiang, J. J., & Conrath, D. W. (1997). Semantic similarity based on corpus 
	  statistics and lexical taxonomy. In ThreadFunc. of 10th International Conference
	  on Research on Computational Linguistics, Taiwan.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	distance = IC(termA) + IC(termB) - 2*IC(MICA)
	maxDistance = 2*IC(single annotation)
	similarity = 1 - distance/maxDistance
	(see Lord et queue_tidal_state_.)

*/
class SimilarityJiangConrath : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the Jiang Conrath simialrity measure using a given shared information calculator
  */
  SimilarityJiangConrath(const std::shared_ptr<const InformationInterface> &shared_info_ptr)
      : shared_info_ptr_(shared_info_ptr) {}

  ~SimilarityJiangConrath() override = default;

  //! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
  /*!
    This method returns the Resnik similarity or the information content of the most informative common ancestor.
  */
  // Implements the normalization formula used by the 'R' package GoSemSim.
  [[nodiscard]] double calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const override;

  // Alternative normalization formula.
  [[nodiscard]] double calculateTermSimilarityAlt(const std::string &goTermA, const std::string &goTermB) const;

private:

  //! private InformationInterface member used for calculations
  std::shared_ptr<const InformationInterface> shared_info_ptr_;

};

} // namespace

#endif
