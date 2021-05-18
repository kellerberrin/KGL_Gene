//
// Created by kellerberrin on 18/5/21.
//

#ifndef KOL_TERMPROBABILITYUNIQUE_H
#define KOL_TERMPROBABILITYUNIQUE_H


#include "kol_GoGraph.h"
#include "kol_AnnotationData.h"
#include "kol_TermProbabilityInterface.h"



namespace kellerberrin::ontology {


/*! \class TermProbabilityMap
	\brief A class to calculate the probability of a GO term.

	This class provides a map that returns the probability of GO term. This
	  class is used by Information Content methods to determine the prior probability of
	  a term give an instance of AnnotationData.
*/


class TermProbabilityUnique : public TermProbabilityInterface {
public:

  //! A default constructor
  /*!
    Default constructor should not be used.
  */
  TermProbabilityUnique() = delete;
  ~TermProbabilityUnique() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized constructor is allowed to ensure these objects are
      created with valid parameters.
  */
  TermProbabilityUnique( const std::shared_ptr<const GoGraph> &graph,
                         const std::shared_ptr<const AnnotationData> &annotation_data) {

    // Set the default minimum probability policy for normalization
    is_single_anno_min_ = false;
    calcProbability(graph, annotation_data);

  }

  //! Accessor for probabilities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] const TermProbOntMap &getValues() const override { return probability_map_; }

  //! Method to test if the id exists in the map
  /*!
    Return true the id is found, false if not
  */
  [[nodiscard]] bool hasTerm(const std::string &testTerm) const override { return probability_map_.contains(testTerm); }



  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] double getValue(const std::string &term_id) const override;

  //-------------------------------------------------------------------

  //! Get the specific minimum probability for BIOLOGICAL_PROCESS
  /*!
    This function returns the minimum probablity for the bp ontology
  */
  [[nodiscard]] double getMinBP() const override {

    return is_single_anno_min_ ? bp_normalization_min_1anno_ : bp_normalization_min_min_anno_;

  }


  //! Get the specific minimum probability for MOLECULAR_FUNCTION
  /*!
    This function returns the minimum probablity for the mf ontology
  */
  [[nodiscard]] double getMinMF() const override {

    return is_single_anno_min_ ? mf_normalization_min_1anno_ : mf_normalization_min_min_anno_;

  }


  //! Get the specific minimum probability for CELLULAR_COMPONENT
  /*!
    This function returns the minimum probablity for the cc ontology
  */
  [[nodiscard]] double getMinCC() const override {

    return is_single_anno_min_ ? cc_normalization_min_1anno_ : cc_normalization_min_min_anno_;

  }


  //! Public method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */

  [[nodiscard]] double getMICAinfo(const OntologySetType<std::string> &ancestorsA,
                                   const OntologySetType<std::string> &ancestorsB) const override;



  [[nodiscard]] TermProbOntMap& probabilityMap() { return probability_map_; }
  // Reset the bad value in information content classes
  void setBadValue(double value) { bad_prob_value_ = value; }

private:

  double bad_prob_value_{1.0};

  //! A private map that returns the index of a term.
  /*!
    This map takes string term ids and returns the index for annotation count access.
  */

  TermProbOntMap probability_map_;

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


  [[nodiscard]] double getRootCount(const std::string& root_id) const;

  void calcProbability(const std::shared_ptr<const GoGraph> &graph,
                       const std::shared_ptr<const AnnotationData> &annotation_data);

  //! Private method for calculating the most informative common ancestor value
  /*!
    This method searches the sets to determine the most informative ancestor.
  */
  [[nodiscard]] double getEfficientMICA(const OntologySetType<std::string> &smaller_set,
                                        const OntologySetType<std::string> &larger_set) const;

  //! Return a default value for a term that does not exist.
  /*!
    A value to return if the term is not found (does not exist in the map).
    Returns probability 1 or certanty. This may not be the ideal behavior.
  */
  [[nodiscard]] double badIdValue() const { return bad_prob_value_; }

};


}  // namespace


#endif //KOL_TERMPROBABILITYUNIQUE_H
