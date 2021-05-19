//
// Created by kellerberrin on 18/5/21.
//

#ifndef KOL_TERMPROBABILITYINTERFACE_H
#define KOL_TERMPROBABILITYINTERFACE_H


#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"



namespace kellerberrin::ontology {


/*! \class TermProbabilityMap
	\brief A class to calculate the probability of a GO term.

	This class provides a map that returns the probability of GO term. This
	  class is used by Information Content methods to determine the prior probability of
	  a term give an instance of AnnotationData.
*/

using TermProbOntMap = OntologyMapType<std::string, std::pair<double, GO::Ontology>>;

class TermProbabilityInterface {
public:

  TermProbabilityInterface() = default;
  virtual ~TermProbabilityInterface() = default;

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] virtual const TermProbOntMap &getValues() const = 0;


  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] virtual double getValue(const std::string &term_id) const = 0;

  //-------------------------------------------------------------------


private:



};


}  // namespace

#endif //KOL_TERMPROBABILITYINTERFACE_H
