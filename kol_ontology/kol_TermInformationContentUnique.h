//
// Created by kellerberrin on 20/5/21.
//

#ifndef KOL_TERMINFORMATIONCONTENTUNIQUE_H
#define KOL_TERMINFORMATIONCONTENTUNIQUE_H


#include "kol_TermInformationContentImpl.h"


namespace kellerberrin::ontology {


/*! \class TermInformationContentMap
	\brief A class to calculate the information content of a GO term.

	This class provides a map that returns the information content of a GO term. This
	  class is used by Information Content methods.

*/
class TermInformationContentUnique : public TermInformationContentImpl {
public:
  //! A default constructor
  /*!
    This constructor creates an empty IC map. Should not be used.
  */
  TermInformationContentUnique() = delete;
  ~TermInformationContentUnique() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized construtor is allowed to ensure these objects are
      created with valid parameters.
      This constructor relies on the TermProbabilityMap.
  */
  TermInformationContentUnique( const std::shared_ptr<const GoGraph> &graph,
                                const std::shared_ptr<const AnnotationData> &annoData) {

    calcProbabilityMap(graph, annoData);
    convertProbtoIC();

  }


private:

  void calcProbabilityMap(const std::shared_ptr<const GoGraph> &graph,
                          const std::shared_ptr<const AnnotationData> &annoData);

};



} // namespace


#endif //KOL_TERMINFORMATIONCONTENTUNIQUE_H
