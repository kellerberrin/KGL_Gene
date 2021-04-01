/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef MODULAR_JIANG_CONRATH
#define MODULAR_JIANG_CONRATH

#include <TermSimilarityInterface.h>
#include <SharedInformationInterface.h>

/*! \class ModularJiangConrath
	\brief A class to calculate Jiang Conrath similarity between 2 terms

	This class calculates Jiang Conrath similarity.
	
	Jiang, J. J., & Conrath, D. W. (1997). Semantic similarity based on corpus 
	  statistics and lexical taxonomy. In Proc. of 10th International Conference
	  on Research on Computational Linguistics, Taiwan.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	distance = IC(termA) + IC(termB) - 2*IC(MICA)
	maxDistance = 2*IC(single annotation)
	similarity = 1 - distance/maxDistance
	(see Lord et al.)

*/
class ModularJiangConrath: public TermSimilarityInterface{

public:
	
	//! A constructor
	/*!
		Creates the Jiang Conrath simialrity measure using a given shared infromation calculator
	*/
	ModularJiangConrath(std::shared_ptr<const SharedInformationInterface> sharedInformationCalculator)
	: _siCalculator(std::move(sharedInformationCalculator)) {}
  ~ModularJiangConrath() override = default;

	//! A method for calculating term-to-term similarity for GO terms using Lin similarity
	/*!
		This method returns the Resnik similarity or the information content of the most informative common ancestor.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

    if (goTermA == goTermB) {

      return 1.0;

    }

		if (not _siCalculator->hasTerm(goTermA) or not _siCalculator->hasTerm(goTermB)){

			return 0.0;

		}

		if (not _siCalculator->isSameOntology(goTermA, goTermB)){

			return 0.0;

		}

		double sharedIC = _siCalculator->sharedInformation(goTermA,goTermB);
		double termA_IC = _siCalculator->sharedInformation(goTermA);
		double termB_IC = _siCalculator->sharedInformation(goTermB);
		double maxIC = _siCalculator->maxInformationContent(goTermA);

		double dist = termA_IC + termB_IC - (2.0 * sharedIC);

		return 1.0 - (dist / (2.0 * maxIC));

	}

	//! A method for calculating term-to-term similarity for GO terms using normalized Lin similarity
	/*!
		This method returns the Lin similarity. Lin similarity is already normalized
	*/
	[[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

		return calculateTermSimilarity(goTermA,goTermB);

	}

private:

	//! private SharedInformationInterface member used for calculations
	std::shared_ptr<const SharedInformationInterface> _siCalculator;

};

#endif
