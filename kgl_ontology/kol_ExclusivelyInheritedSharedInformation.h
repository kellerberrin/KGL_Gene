/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef EXCLUSIVELY_INHERITED_SHARED_INFORMATION
#define EXCLUSIVELY_INHERITED_SHARED_INFORMATION

#include "kol_TermInformationContentMap.h"
#include "kol_SharedInformationInterface.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"
#include "kol_GoGraph.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kellerberrin::ontology {

/*! \class ExclusivelyInheritedSharedInformation
	\brief A class to calculate shared infromation in linear time after Zhang and Lai.

	Shu-Bo Zhang and Jian-Huang Lai. Semantic Similarity measurement between gene ontology
	  terms based on exclusively inherited shared informaiton. Gene 558 (2015) 108-117.

*/
class ExclusivelyInheritedSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMSharedInformation class
  */
  ExclusivelyInheritedSharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                                         const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~ExclusivelyInheritedSharedInformation() override = default;

  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1,
                                                                        const std::string &termC2) const;

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
  [[nodiscard]] bool hasTerm(const std::string &term) const override { return ic_map_ptr_->hasTerm(term); }

  //! An interface method for determining if the two terms are of like ontologies.
  /*!
    Determine if two terms are of the same ontology.
  */
  [[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

    return graph_ptr_->getTermOntology(termA) == graph_ptr_->getTermOntology(termB);

  }

private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationContentMap> ic_map_ptr_;

};

} // namespace

#endif
