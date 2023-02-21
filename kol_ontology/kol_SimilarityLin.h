/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_MODULAR_LIN
#define KOL_MODULAR_LIN


#include "kol_SimilarityInterface.h"
#include "kol_InformationInterface.h"

#include <memory>

namespace kellerberrin::ontology {

/*! \class SimilarityLin
	\brief A class to calculate Lin similarity between 2 terms

	This class calculates Lin similarity.
	
	Lin, D. (1998) An information theoretic definition of similarity. In: ThreadFunc. of the
	  15th Inernational Conference on Machine Learning. San Franscisco, CA:
	  Morgan Kaufman. pp 296-304

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	2 * IC(MICA) / ( IC(termA) + IC(termB) )

*/
class SimilarityLin : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the linSimPtr calculator with a particular shared information calculator
  */
  SimilarityLin(const std::shared_ptr<const InformationInterface> &shared_info_ptr) : shared_info_ptr_(shared_info_ptr) {}
  ~SimilarityLin() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Lin similarity
  /*!
    This method returns the Resnik similarity or the information content of the most informative common ancestor.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const override;


private:

  //! private InformationInterface member used for calculations
  std::shared_ptr<const InformationInterface> shared_info_ptr_;

};


} // namespace

#endif
