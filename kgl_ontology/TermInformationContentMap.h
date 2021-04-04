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

	  if (not checkIntegrity()) {

      throw std::runtime_error("TermProbabilityMap::TermProbabilityMap::failed integrity check.");

	  }
		//convert term probability to information content
    auto index_map = indexToName();
    size_t zero_count{0};
    size_t missing_non_zero{0};
    double prob_sum{0.0};
    double missing_sum{0.0};
		for(size_t i = 0; i < probabilities().size(); ++i){

      prob_sum += probabilities().at(i);
      auto result = index_map.find(i);
      if (result == index_map.end()) {

        throw std::runtime_error("TermInformationContentMap::TermInformationContentMap; index not found");

      }
      auto const& [index, name] = *result;
      if (not annoData->hasGoTerm(name) and probabilities().at(i) != 0) {

        ++missing_non_zero;
        missing_sum += probabilities().at(i);

      }

		  if (probabilities().at(i) == 0.0) {

        ++zero_count;
        probabilities().at(i) = badIdValue();

      } else {

        probabilities().at(i) = -1.0 * std::log(probabilities().at(i));

		  }


		}//end for, each probability value

    std::stringstream ss;
    ss << "TermInformationContentMap::TermInformationContentMap; prob vector size: " << probabilities().size()
       << " zero terms: " << zero_count <<  " missing non-zero: "<< missing_non_zero
       << " missing sum prob: " << missing_sum << " prob sum: " << prob_sum << '\n';
    std::cout << ss.str();

	}

	//! Return a default value for a term that does not exist.
	/*!
	A value to return if the term is not found (does not exist in the map).
	Returns information content 0. This may not be the ideal behavior.
	*/
	double badIdValue() const { return 0.0; }


};
#endif

