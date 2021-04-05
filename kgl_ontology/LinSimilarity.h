/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef LIN_SIMILARITY
#define LIN_SIMILARITY

#include <TermSimilarityInterface.h>
#include <TermInformationContentMap.h>
#include <GoGraph.h>

/*! \class LinSimilarity
	\brief A class to calculate Lin similarity between 2 terms

	This class calculates Lin similarity.
	
	Lin, D. (1998) An information theoretic definition of similarity. In: Proc. of the
	  15th Inernational Conference on Machine Learning. San Franscisco, CA:
	  Morgan Kaufman. pp 296-304

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	2 * IC(MICA) / ( IC(termA) + IC(termB) )

*/
class LinSimilarity : public TermSimilarityInterface{

public:
	
	//! A constructor
	/*!
		Creates the LinSimilarity calculator
	*/
	LinSimilarity(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	:	_goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~LinSimilarity() override = default;

	//! A method for calculating term-to-term similarity for GO terms using Lin similarity
	/*!
		This method returns the Lin similarity.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

    if (goTermA == goTermB) {

      return 1.0;

    }

		//if the terms do not exit return 0.0 similarity
		if (not _icMap->hasTerm(goTermA) || not _icMap->hasTerm(goTermB)) {

			return 0.0;

		}
		//if not from same ontology, return 0;
		if(_goGraph->getTermOntology(goTermA) != _goGraph->getTermOntology(goTermB)){

			return 0.0;

		}

		//create 2 sets
		OntologySetType<std::string> ancestorsA = _goGraph->getAncestorTerms(goTermA);
		ancestorsA.insert(goTermA);
		OntologySetType<std::string> ancestorsB = _goGraph->getAncestorTerms(goTermB);
		ancestorsB.insert(goTermB);

		//if either set is empty, return 0
		if(ancestorsA.empty() or ancestorsB.empty()) {

			return 0.0;

		}

    double mica_info = _icMap->getMICAinfo(ancestorsA,ancestorsB);
    //return the normalized information content similarity of Lin
    return (2.0 * mica_info) /(_icMap->getValue(goTermA) + _icMap->getValue(goTermB));


	}

	//! A method for calculating term-to-term similarity for GO terms using Normalized Lin similarity
	/*!
		This method returns the Lin similarity scaled between 0 and 1 [0,1] inclusive
	*/
	[[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//Lin's method is already normalized
		return calculateTermSimilarity(goTermA,goTermB);

	}

private:

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;

};
#endif
