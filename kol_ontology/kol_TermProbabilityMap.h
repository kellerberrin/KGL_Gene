/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_TERM_PROBABILITY_MAP
#define KGL_TERM_PROBABILITY_MAP


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


class TermProbabilityMap : public TermProbabilityInterface {
public:

  //! A default constructor
  /*!
    Default constructor should not be used.
  */
  TermProbabilityMap() = delete;
  ~TermProbabilityMap() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized constructor is allowed to ensure these objects are
      created with valid parameters.
  */
  TermProbabilityMap(const std::shared_ptr<const GoGraph> &graph,
                     const std::shared_ptr<const AnnotationData> &annotation_data) {

    // Set the default minimum probability policy for normalization
    calcProbabilityMap(graph, annotation_data);

  }

  //! Accessor for probabilities vector
  /*!
    Get the vector of values
  */

  [[nodiscard]] const TermProbOntMap &getValues() const override { return probability_map_; }


  //! Mapping function to return the value mapped by key
  /*!
  Get the value mapped by the given key. A specified function for the [] operator
  */

  [[nodiscard]] double getValue(const std::string &term_id) const override;

  //-------------------------------------------------------------------


  // Only use in the information content object.
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

  void calcProbabilityMap( const std::shared_ptr<const GoGraph> &graph,
                           const std::shared_ptr<const AnnotationData> &annotation_data);


  //! Return a default value for a term that does not exist.
  /*!
    A value to return if the term is not found (does not exist in the map).
    Returns probability 1 or certanty. This may not be the ideal behavior.
  */
  [[nodiscard]] double badIdValue() const { return bad_prob_value_; }

  // Return the number of annotations at the root of each GO ontology.
  [[nodiscard]] double getRootCount(const std::string& root_term) const;

};

}  // namespace

#endif
