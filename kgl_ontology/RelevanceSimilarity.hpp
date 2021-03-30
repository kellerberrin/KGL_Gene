/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef RELEVANCE_SIMILARITY
#define RELEVANCE_SIMILARITY

#include <TermSimilarityInterface.hpp>
#include <TermInformationContentMap.hpp>
#include <GoGraph.hpp>

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
class RelevanceSimilarity: public TermSimilarityInterface {

public:
	
	//! A constructor
	/*!
		Creates the default(empty) StandardRelationshipPolicy
	*/
	RelevanceSimilarity(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	:	_goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~RelevanceSimilarity() override = default;

	//! A method for calculating term-to-term similarity for GO terms using Relevance similarity
	/*!
		This method returns the Relevance similarity.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//if the terms do not exit return 0.0 similarity
		if (not _icMap->hasTerm(goTermA) or not _icMap->hasTerm(goTermB)){
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

		//get the MICA
		double mica_info = _icMap->getMICAinfo(ancestorsA, ancestorsB);
		double complement_prob_mica = 1.0 - std::exp(-1.0 * mica_info);

		//return the normalized information content similarity of Relevance
		return (2.0 * mica_info) / ((_icMap->getValue(goTermA) + _icMap->getValue(goTermB)) * complement_prob_mica);

	}

	//! A method for calculating term-to-term similarity for GO terms using Normalized Relevance similarity
	/*!
		This method returns the Relevance similarity scaled between 0 and 1 [0,1] inclusive
	*/
	[[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//Relevance's method is already normalized
		return calculateTermSimilarity(goTermA,goTermB);

	}


private:

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;

};
#endif