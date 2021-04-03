/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef TERM_INFORMATION_CONTENT_MAP
#define TERM_INFORMATION_CONTENT_MAP

#include <cmath>

#include <GoGraph.h>
#include <AnnotationData.h>
#include <TermProbabilityMap.h>


/*! \class TermInformationContentMap
	\brief A class to calculate the information content of a GO term.

	This class provides a map that returns the information content of a GO term. This
	  class is used by Information Content methods.

*/
class TermInformationContentMap: public TermProbabilityMap{
public:
  //! A default constructor
  /*!
    This constructor creates an empty IC map. Should not be used.
  */
  TermInformationContentMap() = delete;
  TermInformationContentMap(const TermInformationContentMap&) = default;
  ~TermInformationContentMap() override = default;
  //! A parameterized constructor
  /*!
    This constructor takes pointers to GoGraph and AnnotationData objects.
      Only the parameterized construtor is allowed to ensure these objects are
      created with valid parameters.
      This constructor relies on the TermProbabilityMap.
  */
	TermInformationContentMap(const std::shared_ptr<const GoGraph>& graph, const std::shared_ptr<const AnnotationData>& annoData)
	: TermProbabilityMap(graph,annoData)
	{
		//convert term probability to information content
		for(size_t i = 0; i < probabilities().size(); ++i){

		  if (probabilities().at(i) == 0.0) {

		    size_t zero_count{0};
		    double prob_sum{0.0};
		    for (auto const& prob : probabilities()) {

		      if (prob == 0) ++zero_count;
		      prob_sum += prob;

		    }

        std::stringstream ss;
        ss << "TermInformationContentMap::TermInformationContentMap; prob vector size: " << probabilities().size()
           << " zero terms: " << zero_count <<  " at index: " << i << " prob sum: " << prob_sum;
		    throw std::runtime_error(ss.str());

		  }

			probabilities().at(i) = -1.0 * std::log(probabilities().at(i));

		}//end for, each probability value

	}

	//! Return a default value for a term that does not exist.
	/*!
	A value to return if the term is not found (does not exist in the map).
	Returns information content 0. This may not be the ideal behavior.
	*/
	double badIdValue() const { return 0.0; }


};
#endif

