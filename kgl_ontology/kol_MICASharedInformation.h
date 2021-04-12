/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/

#ifndef KGL_MICA_SHARED_INFORMATION
#define KGL_MICA_SHARED_INFORMATION

#include "kol_TermInformationContentMap.h"
#include "kol_SharedInformationInterface.h"
#include "kol_GoGraph.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"

#include <boost/accumulators/statistics/max.hpp>

namespace kellerberrin::ontology {


/*! \class MICASharedInformation
	\brief A class to calculate shared information as the most informative common ancestor (MICA)

	This class calculates shared information using the most informative common ancestor (MICA).
	 The MICA is a term that is also known as the minimum subsumer.

	 This shared information method forms the basis of 3 information content measures
	 put forward by Lord el al.

    P. W. Lord, R. D. Stevens, A. Brass, and C. A. Goble, 
	 "Semantic similarity measures as tools for exploring the gene ontology,"
	 Pac Symp Biocomput, pp. 601-12, 2003.

*/
class MICASharedInformation : public SharedInformationInterface {

public:
  //! A constructor
  /*!
    Creates the MICASharedInformation class
  */
  MICASharedInformation(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~MICASharedInformation() override = default;

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
    if (_goGraph->getTermOntology(termA) != _goGraph->getTermOntology(termB)) {

      return 0.0;

    }

    Accumulators::MaxAccumulator myMax;

    OntologySetType<std::string> ancestorsA = _goGraph->getSelfAncestorTerms(termA);
    OntologySetType<std::string> ancestorsB = _goGraph->getSelfAncestorTerms(termB);

    OntologySetType<std::string> sharedAncestors = SetUtilities::setIntersection(ancestorsA, ancestorsB);

    for (auto const &term : sharedAncestors) {

      myMax(_icMap->getValue(term));

    }

    return Accumulators::extractMax(myMax);

  }

  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] double sharedInformation(const std::string &term) const override {
    // return 0 for any terms not in the database

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


    //select the correct ontology normalization factor
    GO::Ontology ontology = _goGraph->getTermOntology(term);
    double maxIC;

    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        maxIC = _icMap->getMinBP();
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        maxIC = _icMap->getMinMF();
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        maxIC = _icMap->getMinCC();
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        maxIC = 0.0;
        break;

    }

    if (maxIC <= 0.0) {

      return 0.0;

    }

    return -1.0 * std::log(maxIC);

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

} // namespace


#endif
