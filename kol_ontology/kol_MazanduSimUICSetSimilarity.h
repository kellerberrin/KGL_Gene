/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef KGL_MAZANDU_SIMUIC_SET_SIMILARITY
#define KGL_MAZANDU_SIMUIC_SET_SIMILARITY

#include "kol_TermSetSimilarityInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {


/*! \class MazanduSimUICSetSimilarity
\brief A class to calculate Mazandu and Mulder's SimUIC similarity between 2 sets of go terms.

A separate measure from their SimDIC.

Mazandu, G. K., & Mulder, N. J. (2014). Information content-based Gene Ontology functional
similarity measures: which one to use for a given biological data type?. PloS one, 9(12), e113859

*/
class MazanduSimUICSetSimilarity : public TermSetSimilarityInterface {

public:
  //! Constructor
  /*!
  Creates the MazanduSimUICSetSimilarity class assigning the GoGraph private memeber.
  */
  MazanduSimUICSetSimilarity(const std::shared_ptr<const GoGraph> &graph_ptr,
                             const std::shared_ptr<const TermInformationInterface> &ic_map_ptr)
      : graph_ptr_(graph_ptr), ic_map_ptr_(ic_map_ptr) {}

  ~MazanduSimUICSetSimilarity() override = default;


  //! A method for calculating term set to term set similarity for GO terms;
  /*!
  This method returns the best match average similarity.
  */
  [[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &row_terms,
                                           const OntologySetType<std::string> &column_terms) const override;


private:

  //! Pointer to the GoGraph object.
  /*!
  A reference to GO graph to be used.
  */
  std::shared_ptr<const GoGraph> graph_ptr_;

  //! The information content map.
  /*!
  An information content map.
  */
  std::shared_ptr<const TermInformationInterface> ic_map_ptr_;

  //! A method for calculating the sum of information content of the terms in a set.
  /*!
    This method calculates the sum of information content of the terms in a set.
  */
  [[nodiscard]] double calcICSum(const OntologySetType<std::string> &terms) const;

};


} // namespace


#endif
