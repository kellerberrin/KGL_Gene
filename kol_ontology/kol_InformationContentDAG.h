/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_TERM_INFORMATION_CONTENT_MAP
#define KGL_TERM_INFORMATION_CONTENT_MAP


#include "kol_InformationContentImpl.h"


namespace kellerberrin::ontology {


/*! \class InformationContentDAG
	\brief A class to calculate the information content of a GO term.

	This class provides a map that returns the information content of a GO term. This
	  class is used by Information Content methods.

*/
class InformationContentDAG : public InformationContentImpl {
public:
  //! A default constructor
  /*!
    This constructor creates an empty IC map. Should not be used.
  */
  InformationContentDAG() = delete;
  ~InformationContentDAG() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized construtor is allowed to ensure these objects are
      created with valid parameters.
  */
  InformationContentDAG(const std::shared_ptr<const GoGraph> &graph,
                        const std::shared_ptr<const AnnotationData> &annoData) : InformationContentImpl(graph) {

    calcProbabilityMap(graph, annoData);
    convertProbtoIC();

  }


private:

  void calcProbabilityMap(const std::shared_ptr<const GoGraph> &graph,
                          const std::shared_ptr<const AnnotationData> &annoData);

};



} // namespace



#endif

