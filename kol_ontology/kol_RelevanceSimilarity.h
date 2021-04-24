/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_RELEVANCE_SIMILARITY
#define KGL_RELEVANCE_SIMILARITY

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
  RelevanceSimilarity( const std::shared_ptr<const GoGraph> &graph_ptr,
                       const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~RelevanceSimilarity() override = default;

  //! A method for calculating term-to-term similarity for GO terms using Relevance similarity
  /*!
    This method returns the Relevance similarity.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override;

  //! A method for calculating term-to-term similarity for GO terms using Normalized Relevance similarity
  /*!
    This method returns the Relevance similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //Relevance's method is already normalized
    return calculateTermSimilarity(goTermA, goTermB);

  }


private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationContentMap> ic_map_ptr_;

};

}   // namespace

#endif