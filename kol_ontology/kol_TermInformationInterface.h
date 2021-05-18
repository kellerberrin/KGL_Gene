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

class TermInformationInterface {
public:

  TermInformationInterface() = default;
  virtual ~TermInformationInterface() = default;

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] virtual const TermProbOntMap &getValues() const = 0;


  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] virtual bool hasTerm(const std::string &testTerm) const = 0;

  //! Return a default value for a term that does not exist.
  /*!
    A value to return if the term is not found (does not exist in the map).
    Returns probability 1 or certanty. This may not be the ideal behavior.
  */
//  [[nodiscard]] virtual double badIdValue() const = 0;

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] virtual double getValue(const std::string &term_id) const = 0;

  //-------------------------------------------------------------------

  //! Get the specific minimum probability for BIOLOGICAL_PROCESS
  /*!
    This function returns the minimum probablity for the bp ontology
  */
  [[nodiscard]] virtual double getMinBP() const = 0;

  //! Get the specific minimum probability for MOLECULAR_FUNCTION
  /*!
    This function returns the minimum probablity for the mf ontology
  */
  [[nodiscard]] virtual double getMinMF() const = 0;

  //! Get the specific minimum probability for CELLULAR_COMPONENT
  /*!
    This function returns the minimum probablity for the cc ontology
  */
  [[nodiscard]] virtual double getMinCC() const = 0;

  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] virtual double getMICAinfo( const OntologySetType<std::string> &ancestorsA,
                                            const OntologySetType<std::string> &ancestorsB) const = 0;



private:



};


}  // namespace


#endif //KGL_KOL_TERMINFORMATIONINTERFACE_H
