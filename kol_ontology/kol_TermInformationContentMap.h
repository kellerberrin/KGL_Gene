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
    Return true the ids are found and the same ontology, false if not
  */
  [[nodiscard]] bool validateTerms(const std::string &id_termA, const std::string &id_termB) const override;

  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] bool hasTerm(const std::string &id_term) const override { return getValues().contains(id_term); }

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] double getValue(const std::string &term_id) const override { return probability_map_.getValue(term_id); }

  //-------------------------------------------------------------------

  //! Get the maximum information content for an ontology class
  /*!
    This function returns the the maximum information content for an ontology class
  */

  [[nodiscard]] double getMaxInformation(const std::string term_id) const override;

  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] double getMICAinfo(const std::string& go_termA, const std::string& go_termB, const GoGraph &graph) const override;

private:

  ProbMap probability_map_;
  const static constexpr double BAD_INFO_VALUE_{0.0};
  double max_bp_information{0.0};
  double max_mf_information{0.0};
  double max_cc_information{0.0};

  void convertProbtoIC();
  [[nodiscard]] double getEfficientMICA( const OntologySetType<std::string> &smaller_set,
                                         const OntologySetType<std::string> &larger_set) const;
  [[nodiscard]] double getMaxInformation(GO::Ontology ontology) const;

};


/*!
  This function converts the probability terms of the probability map to information content.
*/
template<class ProbMap>
void InformationContentMap<ProbMap>::convertProbtoIC() {

  for (auto& [term_id, prob_ont_pair] : probability_map_.probabilityMap()) {

    auto& [probability, ontology] = prob_ont_pair;

    if (probability <= 0.0) {

      probability = BAD_INFO_VALUE_;

    } else {

      probability = -1.0 * std::log(probability);
      switch(ontology) {

        case GO::Ontology::BIOLOGICAL_PROCESS:
          max_bp_information = std::max(probability, max_bp_information);
          break;

        case GO::Ontology::MOLECULAR_FUNCTION:
          max_mf_information = std::max(probability, max_mf_information);
          break;

        case GO::Ontology::CELLULAR_COMPONENT:
          max_cc_information = std::max(probability, max_cc_information);
          break;

        default:
        case GO::Ontology::ONTO_ERROR:
          break;

      }

    }

  }

  // The probability map has been converted to information content, so change the missing value to 0.0.
  probability_map_.setBadValue(BAD_INFO_VALUE_);

}

// Find if the information map contains two terms and they have the same ontology.
template<class ProbMap>
bool InformationContentMap<ProbMap>::validateTerms(const std::string &id_termA, const std::string &id_termB) const {

  auto result_A = probability_map_.getValues().find(id_termA);
  if (result_A == probability_map_.getValues().end()) {

    return false;

  }

  auto result_B = probability_map_.getValues().find(id_termB);
  if (result_B == probability_map_.getValues().end()) {

    return false;

  }

  auto const& [termA, val_ont_A] = *result_A;
  auto const& [valA, ontA] = val_ont_A;

  auto const& [termB, val_ont_B] = *result_B;
  auto const& [valB, ontB] = val_ont_B;

  return ontA == ontB;

}


//! Public method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/

template<class ProbMap>
double InformationContentMap<ProbMap>::getMICAinfo(const std::string& go_termA, const std::string& go_termB, const GoGraph &graph) const {

  //create 2 sets
  OntologySetType<std::string> ancestorsA = graph.getSelfAncestorTerms(go_termA);
  OntologySetType<std::string> ancestorsB = graph.getSelfAncestorTerms(go_termB);

  if (ancestorsA.empty() or ancestorsB.empty()) {

    return 0.0;

  }

  // Choose the smaller and larger set for maximum efficiency
  if (ancestorsA.size() < ancestorsB.size()) {

    return getEfficientMICA(ancestorsA, ancestorsB);

  } else {

    return getEfficientMICA(ancestorsB, ancestorsA);

  }

}

//! Private method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/
template<class ProbMap>
double InformationContentMap<ProbMap>::getEfficientMICA( const OntologySetType<std::string> &smaller_set,
                                                         const OntologySetType<std::string> &larger_set) const {

  double max{0.0};
  //loop over shorter list
  for (auto const &term : smaller_set) {

    if (larger_set.find(term) != larger_set.end()) {

      double term_value = getValue(term);
      if (term_value > max) {

        max = term_value;

      }

    }

  }

  return max;

}

template<class ProbMap>
double InformationContentMap<ProbMap>::getMaxInformation(GO::Ontology ontology) const {

  switch(ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      return max_bp_information;

    case GO::Ontology::MOLECULAR_FUNCTION:
      return max_mf_information;

    case GO::Ontology::CELLULAR_COMPONENT:
      return max_cc_information;

    default:
    case GO::Ontology::ONTO_ERROR:
      return 0.0;

  }

}


template<class ProbMap>
double InformationContentMap<ProbMap>::getMaxInformation(const std::string term_id) const {

  auto result = probability_map_.getValues().find(term_id);
  if (result == probability_map_.getValues().end()) {

    return BAD_INFO_VALUE_;

  }

  auto const& [term, val_ont_pair] = *result;
  auto const& [value, ontology] = val_ont_pair;

  return getMaxInformation(ontology);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Uses the DAG information/probability content algorithm.
using TermInformationContentMap = InformationContentMap<TermProbabilityMap>;
// Uses the unique descendent term information/probability content algorithm.
using TermInformationContentUnique = InformationContentMap<TermProbabilityUnique>;


} // namespace



#endif

