/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef ANCESTOR_MEAN_SHARED_INFORMATION
#define ANCESTOR_MEAN_SHARED_INFORMATION

#include <TermInformationContentMap.hpp>
#include <SharedInformationInterface.hpp>
#include <GoGraph.hpp>
#include <SetUtilities.hpp>
#include <Accumulators.hpp>


#include <boost/accumulators/statistics/max.hpp>

/*! \class AncestorMeanSharedInformation
	\brief A class to calculate shared infromation as the average information conent of all common ancestors

	This class calculates shared infromation by averaging the information content of all common ancestors.

	 This shared information method is used a baseline for comparison and may not be meaningful.

*/
class AncestorMeanSharedInformation : public SharedInformationInterface{

public:
	//! A constructor
	/*!
		Creates the AncestorMeanSharedInformation class
	*/
	AncestorMeanSharedInformation(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	: _goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~AncestorMeanSharedInformation() override = default;

	//! A method for calculating the shared infromation between two concepts.
	/*!
		This method returns the shared information between two concepts.
	*/
	[[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override {
		// return 0 for any terms not in the datbase
		if (not _icMap->hasTerm(termA) or not _icMap->hasTerm(termB)) {

			return 0.0;

		}
		// return 0 for terms in different ontologies
		if (not isSameOntology(termA, termB)){

			return 0.0;

		}

		Accumulators::MeanAccumulator meanIC;

    OntologySetType<std::string> ancestorsA = _goGraph->getAncestorTerms(termA);
		ancestorsA.insert(termA);
    OntologySetType<std::string> ancestorsB = _goGraph->getAncestorTerms(termB);
		ancestorsB.insert(termB);

		OntologySetType<std::string> sharedAncestors = SetUtilities::set_intersection(ancestorsA,ancestorsB);

		for(auto const& term : sharedAncestors) {

			meanIC(_icMap->getValue(term));

		}

		return Accumulators::extractMean(meanIC);

	}

	//! An interface method for returning the shared information of a single terms,or information content
	/*!
		This method privdes a mechanism for returing a term's infromation content.
	*/
	[[nodiscard]] double sharedInformation(const std::string &term) const override {
		// return 0 for any terms not in the datbase
		if (not _icMap->hasTerm(term)) {

			return 0.0;

		}

		return _icMap->getValue(term);

	}


	//! An interface method for returning the maximum information content for a term
	/*!
		This method provides the absolute max information content within a corpus for normalization purposes.
	*/
	[[nodiscard]] double maxInformationContent(const std::string &term) const override {

		double maxIC;

		//select the correct ontology normalization factor
		GO::Ontology ontoType = _goGraph->getTermOntology(term);
		if(ontoType == GO::Ontology::BIOLOGICAL_PROCESS){

			maxIC = -std::log(_icMap->getMinBP());

		}else if(ontoType == GO::Ontology::MOLECULAR_FUNCTION){

			maxIC = -std::log(_icMap->getMinMF());

		}else{

			maxIC = -std::log(_icMap->getMinCC());
		}

		return maxIC;

	}

	//! An interface method for determining if a term can be found
	/*!
		Determines if the term can be found in the current map.
	*/
	[[nodiscard]] bool hasTerm(const std::string &term) const override { return _icMap->hasTerm(term); }

	//! An interface method for determining if the two terms are of like ontologies.
	/*!
		Determine if two terms are of the same ontology.
	*/
	[[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

		return _goGraph->getTermOntology(termA) == _goGraph->getTermOntology(termB);

	}

private:

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;


};
#endif
