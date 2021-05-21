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

	This class calculates shared information along a semantic frontier between terms.
*/
class FrontierSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMSharedInformation class
  */
  FrontierSharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                             const std::shared_ptr<const TermInformationInterface> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~FrontierSharedInformation() override = default;


  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content of the frontier ancestors
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override;
  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] double termInformation(const std::string &term) const override;
  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override;

  //! Method to test if the ids exist and have the same ontology in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] bool validateTerms(const std::string &id_termA, const std::string &id_termB) const override {

    return ic_map_ptr_->validateTerms(id_termA, id_termB);

  }

private:

  std::shared_ptr<const GoGraph> graph_ptr_;
  std::shared_ptr<const TermInformationInterface> ic_map_ptr_;

  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const;

};

} // namespace

#endif
