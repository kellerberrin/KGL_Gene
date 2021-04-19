/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_TERM_PROBABILITY_MAP
#define KGL_TERM_PROBABILITY_MAP

#include <cmath>
#include <string>

#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_Accumulators.h"



namespace kellerberrin::ontology {


/*! \class TermProbabilityMap
	\brief A class to calculate the probability of a GO term.

	This class provides a map that returns the probability of GO term. This
	  class is used by Information Content methods to determine the prior probability of
	  a term give an instance of AnnotationData.
*/
class TermProbabilityMap {
public:

  //! A default constructor
  /*!
    Default constructor should not be used.
  */
  TermProbabilityMap() = delete;
  virtual ~TermProbabilityMap() = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized constructor is allowed to ensure these objects are
      created with valid parameters.
  */
  TermProbabilityMap(const std::shared_ptr<const GoGraph> &graph,
                     const std::shared_ptr<const AnnotationData> &annoData);

  //! Accessor for probablities vector
  /*!
    Get the vector of values
  */
  [[nodiscard]] const std::vector<double> &getValues() const { return probabilities_; }

  //! Function to return all the keys in the map
  /*!
    Returns all valid keys in the map.
  */
  [[nodiscard]] std::vector<std::string> getKeys() const;

  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] bool hasTerm(const std::string &testTerm) const { return name_to_index_.find(testTerm) != name_to_index_.end(); }

  //! Return a default value for a term that does not exist.
  /*!
    A value to return if the term is not found (does not exist in the map).
    Returns probability 1 or certanty. This may not be the ideal behavior.
  */
  [[nodiscard]] virtual double badIdValue() const { return BAD_PROB_VALUE_; }

  //! Overloaded [] bracket operator to mimic Map
  /*!
  This defines a bracket operator to access the data inside of the map.
  This is done to mimic the behavior of the map class
  */
  [[nodiscard]] double operator[](const std::string &termId) const {

    return getValue(termId);

  }

  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */
  [[nodiscard]] double getValue(const std::string &termId) const;

  //-------------------------------------------------------------------

  //! Get the specific minimum probability for BIOLOGICAL_PROCESS
  /*!
    This function returns the minimum probablity for the bp ontology
  */
  [[nodiscard]] double getMinBP() const {

    return is_single_anno_min_ ? bp_normalization_min_1anno_ : bp_normalization_min_min_anno_;

  }


  //! Get the specific minimum probability for MOLECULAR_FUNCTION
  /*!
    This function returns the minimum probablity for the mf ontology
  */
  [[nodiscard]] double getMinMF() const {

    return is_single_anno_min_ ? mf_normalization_min_1anno_ : mf_normalization_min_min_anno_;

  }


  //! Get the specific minimum probability for CELLULAR_COMPONENT
  /*!
    This function returns the minimum probablity for the cc ontology
  */
  [[nodiscard]] double getMinCC() const {

    return is_single_anno_min_ ? cc_normalization_min_1anno_ : cc_normalization_min_min_anno_;

  }


  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] double getMICAinfo(const OntologySetType<std::string> &ancestorsA,
                                   const OntologySetType<std::string> &ancestorsB) const;


protected:

  [[nodiscard]] std::vector<double> &probabilities() { return probabilities_; }

private:

  const static constexpr double BAD_PROB_VALUE_{1.0};

  //! A private map that returns the index of a term.
  /*!
    This map takes string term ids and returns the index for annotation count access.
  */
  OntologyMapType<std::string, std::size_t> name_to_index_;

  //! A private list of term probabilities
  /*!
    This vector of doubles holds the prior probability for each term
  */
  std::vector<double> probabilities_;

  //! A flag designating the minimum policy
  /*!
    This flag will be true and return true is single annotation probability is used, false otherwise.
  */
  bool is_single_anno_min_{true};

  //! Normalization factor for calculating normalized similarities Biological Process
  /*!
    Normalization factor representing the minimum probability using a single annotation
      divided by the cumulative annotations.
  */
  double bp_normalization_min_1anno_;

  //! Normalization factor for calculating normalized simialrites for Biological Process
  /*!
    Normalization factor representing the minimum probability using the number of
      annotations of the least probable term devided by the cumulative annotations.
  */
  double bp_normalization_min_min_anno_;

  //! Normalization factor for calculating normalized simialrites Molecular Function
  /*!
    Normalization factor representing the minimum probability using a single annotation
      devided by the cumulative annotations.
  */
  double mf_normalization_min_1anno_;

  //! Normalization factor for calculating normalized simialrites for Molecular Function
  /*!
    Normalization factor representing the minimum probability using the number of
      annotations of the least probable term devided by the cumulative annotations.
  */
  double mf_normalization_min_min_anno_;

  //! Normalization factor for calculating normalized simialrites Cellular Component
  /*!
    Normalization factor representing the minimum probability using a single annotation
      devided by the cumulative annotations.
  */
  double cc_normalization_min_1anno_;

  //! Normalization factor for calculating normalized simialrites for Cellular Component
  /*!
    Normalization factor representing the minimum probability using the number of
      annotations of the least probable term devided by the cumulative annotations.
  */
  double cc_normalization_min_min_anno_;


  //! Private method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */
  [[nodiscard]] double getEfficientMICA(const OntologySetType<std::string> &smaller_set,
                                        const OntologySetType<std::string> &larger_set) const;


};

}  // namespace

#endif
