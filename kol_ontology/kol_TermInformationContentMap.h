/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_TERM_INFORMATION_CONTENT_MAP
#define KGL_TERM_INFORMATION_CONTENT_MAP

#include <cmath>

#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermProbabilityMap.h"
#include "kol_TermProbabilityUnique.h"
#include "kol_TermInformationInterface.h"


namespace kellerberrin::ontology {


/*! \class TermInformationContentMap
	\brief A class to calculate the information content of a GO term.

	This class provides a map that returns the information content of a GO term. This
	  class is used by Information Content methods.

*/
template<class ProbMap>
class InformationContentMap : public TermInformationInterface {
public:
  //! A default constructor
  /*!
    This constructor creates an empty IC map. Should not be used.
  */
  InformationContentMap() = delete;

  InformationContentMap(const InformationContentMap &) = default;

  ~InformationContentMap() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized construtor is allowed to ensure these objects are
      created with valid parameters.
      This constructor relies on the TermProbabilityMap.
  */
  InformationContentMap( const std::shared_ptr<const GoGraph> &graph,
                         const std::shared_ptr<const AnnotationData> &annoData)
                         : probability_map_(graph, annoData) {

    convertProbtoIC();

  }

  //! Return a default value for a term that does not exist.
  /*!
  A value to return if the term is not found (does not exist in the map).
  Returns information content 0. This may not be the ideal behavior.
  */

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] const TermProbOntMap &getValues() const override { return probability_map_.getValues(); }

  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] bool hasTerm(const std::string &id_term) const override { return probability_map_.hasTerm(id_term); }

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] double getValue(const std::string &term_id) const override { return probability_map_.getValue(term_id); }

  //-------------------------------------------------------------------

  //! Get the specific minimum probability for BIOLOGICAL_PROCESS
  /*!
    This function returns the minimum probablity for the bp ontology
  */
  [[nodiscard]] double getMinBP() const override { return probability_map_.getMinBP(); }

  //! Get the specific minimum probability for MOLECULAR_FUNCTION
  /*!
    This function returns the minimum probablity for the mf ontology
  */
  [[nodiscard]] double getMinMF() const override { return probability_map_.getMinMF(); }

  //! Get the specific minimum probability for CELLULAR_COMPONENT
  /*!
    This function returns the minimum probablity for the cc ontology
  */
  [[nodiscard]] double getMinCC() const override { return probability_map_.getMinCC(); }

  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] double getMICAinfo( const OntologySetType<std::string> &ancestorsA,
                                    const OntologySetType<std::string> &ancestorsB) const override {

    return probability_map_.getMICAinfo(ancestorsA, ancestorsB);

  }


private:

  ProbMap probability_map_;
  const static constexpr double BAD_INFO_VALUE_{0.0};

  void convertProbtoIC();

};


//! A parameterized constructor
/*!
  This constructor takes pointers to GoGraph and AnnotationData objects.
    Only the parameterized construtor is allowed to ensure these objects are
    created with valid parameters.
    This constructor relies on the TermProbabilityMap.
*/
template<class ProbMap>
void InformationContentMap<ProbMap>::convertProbtoIC() {

  for (auto& [term_id, prob_ont_pair] : probability_map_.probabilityMap()) {

    auto& [probability, ontology] = prob_ont_pair;

    if (probability <= 0.0) {

      probability = BAD_INFO_VALUE_;

    } else {

      probability = -1.0 * std::log(probability);

    }

  }

  probability_map_.setBadValue(0.0);

}


// Uses the DAG information/probability content algorithm.
using TermInformationContentMap = InformationContentMap<TermProbabilityMap>;
// Uses the unique descendent term information/probability content algorithm.
using TermInformationContentUnique = InformationContentMap<TermProbabilityUnique>;


} // namespace



#endif

