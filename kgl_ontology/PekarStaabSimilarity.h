/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef PEKAR_STAAB_SIMILARITY
#define PEKAR_STAAB_SIMILARITY

#include <TermSimilarityInterface.h>
#include <TermDepthMap.h>
#include <GoGraph.h>

/*! \class PekarStaabSimilarity
	\brief A class to calculate PekarStaab similarity between 2 terms

	This class calculates Pekar Staab similarity.
	
	V. Pekar and S. Staab, "Taxonomy learning: factoring the structure 
	 of a taxonomy into a semantic classification decision," in 
	 Proc. of 19th International Conference on Computational Linguistics. 
	 Morristown NJ USA: Association for Computational Linguistics, pp. 1-7, 2002.

	H. Yu, L. Gao, K. Tu, and Z. Guo, "Broadly predicting specific gene 
	 functions with expression similarity and taxonomy similarity,"
	 Gene, vol. 352, pp. 75-81, Jun 6 2005.
	  
    lowest common ancestor (LCA)
	GraphDist(LCA,root)/(GraphDist(a,LCA)+GraphDist(b,LCA)+GraphDist(LCA,root))

*/
class PekarStaabSimilarity: public TermSimilarityInterface {

public:
	
	//! A constructor
	/*!
		Creates the default(empty) StandardRelationshipPolicy
	*/
	PekarStaabSimilarity(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermDepthMap> icMap)
	:	_goGraph(std::move(goGraph)), _depthMap(std::move(icMap)) {}
  ~PekarStaabSimilarity() override = default;


	//! A method for calculating term-to-term similarity for GO terms using Pekar Staab similarity
	/*!
		This method returns the PekarStaab similarity.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//if the terms do not exit return 0.0 similarity
		if (not _depthMap->hasTerm(goTermA) or not _depthMap->hasTerm(goTermB)) {

			return 0.0;

		}

		//if not from same ontology, return 0;
		if(_goGraph->getTermOntology(goTermA) != _goGraph->getTermOntology(goTermB)){

			return 0.0;

		}

		//Create 2 sets self + ancestors.
		OntologySetType<std::string> ancestorsA = _goGraph->getSelfAncestorTerms(goTermA);
		OntologySetType<std::string> ancestorsB = _goGraph->getSelfAncestorTerms(goTermB);

		//if either set is empty, return 0
		if(ancestorsA.empty() || ancestorsB.empty()) {

			return 0.0;

		}

		std::string lca = _depthMap->getLCA(ancestorsA, ancestorsB);
		TermDepthType lcaDepth = _depthMap->getValue(lca);
    TermDepthType denom = (_depthMap->getValue(goTermA) - lcaDepth) + (_depthMap->getValue(goTermB) - lcaDepth) + lcaDepth;

    if (denom == 0) {

      return 0;

    } else {

      return static_cast<double>(lcaDepth) / static_cast<double>(denom);

    }

	}

	//! A method for calculating term-to-term similarity for GO terms using Normalized Pekar Staab similarity
	/*!
		This method returns the PekarStaab similarity scaled between 0 and 1 [0,1] inclusive
	*/
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const {
		//if the terms do not exit return 0.0 similarity
		if (not _depthMap->hasTerm(goTermA) || not _depthMap->hasTerm(goTermB)) {

			return 0.0;

		}
		//Pekar and Staab's method is already normalized
		return calculateTermSimilarity(goTermA,goTermB);

	}


private:

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermDepthMap> _depthMap;

};
#endif