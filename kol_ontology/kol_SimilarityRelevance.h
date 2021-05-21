//
// Created by kellerberrin on 19/5/21.
//

#ifndef KOL_MODULARRELEVANCE_H
#define KOL_MODULARRELEVANCE_H


#include "kol_TermSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {

//! A class to calculate Relevance similarity between 2 terms
/*! \class RelevanceSimilarity

	This class calculates Relevance similarity.

	A. Schlicker, F. S. Domingues, J. Rahnenfuhrer, and T. Lengauer,
	 "A new measure for functional similarity of gene products based
	 on Gene Ontology," BMC Bioinformatics, vol. 7, p. 302, 2006.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble,
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

	  Basically this is Lin similarity scaled by the
	   complement of the probability of the mica
	2 * IC(MICA) / ( IC(termA) + IC(termB) )*(1-p(Mica))
*/
class RelevanceSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  RelevanceSimilarity(const std::shared_ptr<const SharedInformationInterface> &shared_info_ptr)
      : shared_info_ptr_(shared_info_ptr) {}

  ~RelevanceSimilarity() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Relevance similarity
  /*!
    This method returns the Relevance similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &go_termA, const std::string &go_termB) const override;


private:

  //! private SharedInformationInterface member used for calculations
  std::shared_ptr<const SharedInformationInterface> shared_info_ptr_;

};

}   // namespace

#endif //KOL_MODULARRELEVANCE_H
