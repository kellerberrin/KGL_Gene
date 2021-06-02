/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KOL_GENTLEMAN_SIMUI_SET_SIMILARITY
#define KOL_GENTLEMAN_SIMUI_SET_SIMILARITY

#include "kol_SetSimilarityInterface.h"
#include "kol_SetUtilities.h"
#include "kol_GoGraph.h"


namespace kellerberrin::ontology {


/*! \class GentlemanSimUISetSimilarity
	\brief A class to calculate Gentleman's UI similarity between go terms for 2 sets.

	Gentlman R. Visualizing and Distances Using GO. URL http://www.bioconductor.org/docs/vignettes.html.

*/
class GentlemanSimUISetSimilarity : public SetSimilarityInterface {

public:
  //! Constructor
  /*!
    Creates the GentlemanUISimilarity class assigning the GoGraph private member.
  */
  explicit GentlemanSimUISetSimilarity(const std::shared_ptr<const GoGraph> &graph_ptr) : graph_ptr_(graph_ptr) {}
  ~GentlemanSimUISetSimilarity() override = default;


  //! A method for calculating term set to term set similarity for GO terms;
  /*!
    This method returns the best match average similarity.
  */
  [[nodiscard]] double calculateSimilarity( const OntologySetType<std::string> &row_terms,
                                            const OntologySetType<std::string> &column_terms) const override;

private:

  //! Pointer to the GoGraph object
  /*!
    A reference to GO graph to be used.
  */
  std::shared_ptr<const GoGraph> graph_ptr_;


};

} // namerspace

#endif
