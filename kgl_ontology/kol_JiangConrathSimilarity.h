/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_JIANG_CONRATH_SIMILARITY
#define KGL_JIANG_CONRATH_SIMILARITY

#include "kol_TermSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include <kol_GoGraph.h>

namespace kellerberrin::ontology {


/*! \class JiangConrathSimilarity
	\brief A class to calculate Jiang Conrath similarity between 2 terms

	This class calculates Jiang Conrath similarity.
	
	Jiang, J. J., & Conrath, D. W. (1997). Semantic similarity based on corpus 
	  statistics and lexical taxonomy. In Proc. of 10th International Conference
	  on Research on Computational Linguistics, Taiwan.

	P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.
	  
	distance = IC(termA) + IC(termB) - 2*IC(MICA)
	maxDistance = 2*IC(single annotaiotn)
	similarity = 1 - distance/maxDistance
	(see Lord et al.)

*/
class JiangConrathSimilarity : public TermSimilarityInterface {

public:

  //! A constructor
  /*!
    Creates the default(empty) StandardRelationshipPolicy
  */
  JiangConrathSimilarity(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~JiangConrathSimilarity() override = default;
  //! A method for calculating term-to-term similarity for GO terms using JiangConrath similarity
  /*!
    This method returns the JiangConrath similarity or the information content of the most informative common ancestor.
  */
  [[nodiscard]] double calculateTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //if the terms do not exit return 0.0 similarity
    if (not _icMap->hasTerm(goTermA) or not _icMap->hasTerm(goTermB)) {

      return 0.0;

    }

    //if not from same ontology, return 0;
    if (_goGraph->getTermOntology(goTermA) != _goGraph->getTermOntology(goTermB)) {

      return 0;

    }

    //create 2 sets
    OntologySetType<std::string> ancestorsA = _goGraph->getAncestorTerms(goTermA);
    ancestorsA.insert(goTermA);
    OntologySetType<std::string> ancestorsB = _goGraph->getAncestorTerms(goTermB);
    ancestorsB.insert(goTermB);

    //if either set is empty, return 0
    if (ancestorsA.empty() or ancestorsB.empty()) {

      return 0.0;

    }

    double maxIC;
    //select the correct ontology normalization factor
    GO::Ontology ontoType = _goGraph->getTermOntology(goTermA);
    if (ontoType == GO::Ontology::BIOLOGICAL_PROCESS) {

      maxIC = -std::log(_icMap->getMinBP());

    } else if (ontoType == GO::Ontology::MOLECULAR_FUNCTION) {

      maxIC = -std::log(_icMap->getMinMF());

    } else {

      maxIC = -std::log(_icMap->getMinCC());

    }

    //get the MICA value (zero if no mica term)
    double mica_value = _icMap->getMICAinfo(ancestorsA, ancestorsB);

    double dist = _icMap->getValue(goTermA) + _icMap->getValue(goTermB) - (2 * mica_value);

    return 1 - (dist / (2.0 * maxIC));
  }

  //! A method for calculating term-to-term similarity for GO terms using Normalized JiangConrath similarity
  /*!
    This method returns the JiangConrath similarity scaled between 0 and 1 [0,1] inclusive
  */
  [[nodiscard]] double calculateNormalizedTermSimilarity(const std::string &goTermA, const std::string &goTermB) const override {
    //if the terms do not exit return 0.0 similarity
    if (not _icMap->hasTerm(goTermA) || not _icMap->hasTerm(goTermB)) {

      return 0.0;

    }

    //JiangConrath's method is already normalized
    return calculateTermSimilarity(goTermA, goTermB);

  }


private:

  std::shared_ptr<const GoGraph> _goGraph;
  std::shared_ptr<const TermInformationContentMap> _icMap;


};


} // namespace


#endif
