//
// Created by kellerberrin on 18/5/21.
//

#ifndef KGL_KOL_TERMINFORMATIONINTERFACE_H
#define KGL_KOL_TERMINFORMATIONINTERFACE_H


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

class TermInformationInterface {
public:

  TermInformationInterface() = default;
  virtual ~TermInformationInterface() = default;

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] virtual const TermProbOntMap &getValues() const = 0;

  //! Method to test if the ids exist and have the same ontology in the map
  /*!
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] virtual bool validateTerms(const std::string &id_termA, const std::string &id_termB) const = 0;

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] virtual double getValue(const std::string &term_id) const = 0;

  //-------------------------------------------------------------------

  //! Get the maximum information content for an ontology
  /*!
    This function returns the the maximum information content for an ontology
  */

  [[nodiscard]] virtual double getMaxInformation(const std::string& term_id) const = 0;

  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] virtual double getMICAinfo(const std::string& go_termA, const std::string& go_termB, const GoGraph &graph) const = 0;


private:



};


}  // namespace


#endif //KGL_KOL_TERMINFORMATIONINTERFACE_H
