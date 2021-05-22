/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_MODULAR_RESNIK
#define KGL_MODULAR_RESNIK

#include "kol_SimilarityInterface.h"
#include "kol_InformationInterface.h"

namespace kellerberrin::ontology {

/*! \class SimilarityResnik
	\brief A class to calculate resnik similarity between 2 terms using a shared information interface

	This class calculates Resnik similarity.
	Philip Resnik (1995). "Using information content to evaluate semantic similarity in a taxonomy". 
	  In Chris S. Mellish (Ed.). Proceedings of the 14th international joint conference on Artificial intelligence (IJCAI'95)

    P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

	maximun information content of all shared ancestors
	IC(MICA)

*/
class SimilarityResnik : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  SimilarityResnik(const std::shared_ptr<const InformationInterface> &shared_info_ptr) : shared_info_ptr_(shared_info_ptr) {}
  ~SimilarityResnik() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Resnik similarity
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
