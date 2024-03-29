//
// Created by kellerberrin on 19/5/21.
//

#ifndef KOL_MODULARRELEVANCE_H
#define KOL_MODULARRELEVANCE_H


#include "kol_SimilarityInterface.h"
#include "kol_InformationContentDAG.h"


namespace kellerberrin::ontology {

//! A class to calculate Relevance similarity between 2 terms
/*! \class SimilarityRelevance

	This class calculates Relevance similarity.

	A. Schlicker, F. S. Domingues, J. Rahnenfuhrer, and T. Lengauer,
	 "A new measure for functional similarity of gene products based
	 on Gene Ontology," BMC Bioinformatics, vol. 7, p. 302, 2006.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble,
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

	  Basically this is Lin similarity scaled by the
	   complement of the probfailure of the mica
	2 * IC(MICA) / ( IC(termA) + IC(termB) )*(1-p(Mica))
*/
class SimilarityRelevance : public SimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) PolicyStandardRelationship
  */
  SimilarityRelevance(const std::shared_ptr<const InformationInterface> &shared_info_ptr)
      : shared_info_ptr_(shared_info_ptr) {}

  ~SimilarityRelevance() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Relevance similarity
  /*!
    This method returns the Relevance similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const override;


private:

  //! private InformationInterface member used for calculations
  std::shared_ptr<const InformationInterface> shared_info_ptr_;

};

}   // namespace

#endif //KOL_MODULARRELEVANCE_H
