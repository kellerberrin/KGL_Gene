/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef MODULAR_LIN
#define MODULAR_LIN

#include <TermSimilarityInterface.h>
#include <SharedInformationInterface.h>

/*! \class ModularLin
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
class ModularLin: public TermSimilarityInterface {

public:
	
	//! A constructor
	/*!
		Creates the LinSimilarity calculator with a particular shared information calculator
	*/
	ModularLin(std::shared_ptr<const SharedInformationInterface> sharedInformationCalculator)
	:	_siCalculator(std::move(sharedInformationCalculator)) {}
  ~ModularLin() override = default;

	//! A method for calculating term-to-term similarity for GO terms using Lin similarity
	/*!
		This method returns the Resnik similarity or the information content of the most informative common ancestor.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

    if(goTermA == goTermB) {

      return 1.0;

    }

		if (not _siCalculator->hasTerm(goTermA) or not _siCalculator->hasTerm(goTermB)){

			return 0.0;

		}

		if (not _siCalculator->isSameOntology(goTermA, goTermB)) {

			return 0.0;

		}

		double sharedIC = _siCalculator->sharedInformation(goTermA,goTermB);
		double termA_IC = _siCalculator->sharedInformation(goTermA);
		double termB_IC = _siCalculator->sharedInformation(goTermB);

		return (2.0 * sharedIC) /(termA_IC + termB_IC);

	}


	//! A method for calculating term-to-term similarity for GO terms using normalized Lin similarity
	/*!
		This method returns the Lin similarity. Lin similarity is already normalized
	*/
	[[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {

		return calculateTermSimilarity(goTermA, goTermB);

	}


private:

	//! private SharedInformationInterface member used for calculations
	std::shared_ptr<const SharedInformationInterface> _siCalculator;

};
#endif
