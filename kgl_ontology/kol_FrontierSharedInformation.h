/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_FRONTIER_SHARED_INFORMATION
#define KGL_FRONTIER_SHARED_INFORMATION

#include "kol_TermInformationContentMap.h"
#include "kol_SharedInformationInterface.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"
#include "kol_GoGraph.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kellerberrin::ontology {


/*! \class FrontierSharedInformation
	\brief A class to calculate shared information across disjoint common ancestors in linear time.

	This class calculates shared infromation along a semantic frontier between terms.
*/
class FrontierSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMSharedInformation class
  */
  FrontierSharedInformation(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~FrontierSharedInformation() override = default;


  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const;
  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content of the frontier ancestors
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;
  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] double sharedInformation(const std::string &term) const override;
  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override;
  //! An interface method for determining if a term can be found
  /*!
    Determines if the term can be found in the current map.
  */
  [[nodiscard]] bool hasTerm(const std::string &term) const override {

    return _icMap->hasTerm(term);

  }

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
