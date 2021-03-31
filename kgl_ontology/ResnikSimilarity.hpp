#ifndef RESNIK_SIMILARITY
#define RESNIK_SIMILARITY

#include <cmath>
#include <string>

#include <TermSimilarityInterface.hpp>
#include <TermInformationContentMap.hpp>
#include <GoGraph.hpp>

#include <boost/unordered_set.hpp>

/*! \class ResnikSimilarity
	\brief A class to calculate resnik similarity between 2 terms

	This class calculates Resnik similarity.
	Philip Resnik (1995). "Using information content to evaluate semantic similarity in a taxonomy". 
	  In Chris S. Mellish (Ed.). Proceedings of the 14th international joint conference on Artificial intelligence (IJCAI'95)

    P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

	maximun information content of all shared ancestors
	IC(MICA)

*/
class ResnikSimilarity: public TermSimilarityInterface{

public:
	
	//! A constructor
	/*!
		Creates the default(empty) StandardRelationshipPolicy
	*/
	ResnikSimilarity(std::shared_ptr<const GoGraph> goGraph, std::shared_ptr<const TermInformationContentMap> icMap)
	:	_goGraph(std::move(goGraph)), _icMap(std::move(icMap)) {}
  ~ResnikSimilarity() override = default;

	//! A method for calculating term-to-term similarity for GO terms using Resnik similarity
	/*!
		This method returns the Resnik similarity or the information content of the most informative common ancestor.
	*/
	[[nodiscard]] double calculateTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//if the terms do not exit return 0.0 similarity
		if(not _icMap->hasTerm(goTermA) or not _icMap->hasTerm(goTermB)) {

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

		//return the information content of the mica
		return _icMap->getMICAinfo(ancestorsA, ancestorsA);

	}

	//! A method for calculating term-to-term similarity for GO terms using Normalized Resnik similarity
	/*!
		This method returns the Resnik similarity divided by the maximum possible similarity
	*/
	[[nodiscard]] double calculateNormalizedTermSimilarity(const std::string& goTermA, const std::string& goTermB) const override {
		//if the terms do not exit return 0.0 similarity
		if (not _icMap->hasTerm(goTermA) or not _icMap->hasTerm(goTermB)) {

			return 0.0;

		}

		//call base similarity
		double resnik = calculateTermSimilarity(goTermA,goTermB);

		double maxIC;
		//select the correct ontology normalization factor
		GO::Ontology ontoType = _goGraph->getTermOntology(goTermA);
		if(ontoType == GO::Ontology::BIOLOGICAL_PROCESS){

			maxIC = -std::log(_icMap->getMinBP());

		}else if(ontoType == GO::Ontology::MOLECULAR_FUNCTION){

			maxIC = -std::log(_icMap->getMinMF());

		}else{

			maxIC = -std::log(_icMap->getMinCC());

		}

		return resnik / maxIC;
	}

private:

	std::shared_ptr<const GoGraph> _goGraph;
	std::shared_ptr<const TermInformationContentMap> _icMap;

};
#endif
