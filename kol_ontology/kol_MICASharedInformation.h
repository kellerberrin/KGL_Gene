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
  MICASharedInformation( const std::shared_ptr<const GoGraph> &graph_ptr,
                         const std::shared_ptr<const TermInformationContentMap> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~MICASharedInformation() override = default;

  //! A method for calculating the shared infromation between two concepts.
  /*!
    This method returns the shared information between two concepts.
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
