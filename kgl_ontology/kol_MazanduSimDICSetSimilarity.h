/*=============================================================================
Copyright (c) 2016 Paul W. Bible

Distributed under the Boost Software License, Version 1.0. (See accompanying
file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
==============================================================================*/
#ifndef MAZANDU_SIMDIC_SET_SIMILARITY
#define MAZANDU_SIMDIC_SET_SIMILARITY

#include <kol_TermSetSimilarityInterface.h>
#include <kol_TermInformationContentMap.h>
#include <kol_SetUtilities.h>
#include <kol_GoGraph.h>


namespace kellerberrin::ontology {

/*! \class MazanduSimDICSetSimilarity
\brief A class to calculate Mazandu and Mulder's SimDIC similarity between 2 sets of go terms.

Mazandu, G. K., & Mulder, N. J. (2014). Information content-based Gene Ontology functional
similarity measures: which one to use for a given biological data type?. PloS one, 9(12), e113859

*/
class MazanduSimDICSetSimilarity : public TermSetSimilarityInterface {

public:
  //! Constructor
  /*!
  Creates the MazanduSimGICSetSimilarity class assigning the GoGraph private member.
  */
  MazanduSimDICSetSimilarity(const std::shared_ptr<const GoGraph> &graph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _graph(graph), _icMap(icMap) {}

  ~MazanduSimDICSetSimilarity() override = default;

  //! A method for calculating term set to term set similarity for GO terms;
  /*!
  This method returns the best match average similarity.
  */
  [[nodiscard]] double calculateSimilarity(const OntologySetType<std::string> &termsA, const OntologySetType<std::string> &termsB) const override {

    // Get the induced set of terms for each set
    OntologySetType<std::string> inducedTermSetA = _graph->getExtendedTermSet(termsA);
    OntologySetType<std::string> inducedTermSetB = _graph->getExtendedTermSet(termsB);
    // Calculate union and intersection
    OntologySetType<std::string> intersection_set = SetUtilities::set_intersection(inducedTermSetA, inducedTermSetB);

    double intersection_sum = calcICSum(intersection_set);
    double setA_sum = calcICSum(inducedTermSetA);
    double setB_sum = calcICSum(inducedTermSetB);

    //if the union is 0, return 0. No division by 0.
    if (setA_sum + setB_sum == 0.0) {

      return 0.0;

    } else {

      return 2 * intersection_sum / (setA_sum + setB_sum);

    }

  }

private:

  //! Pointer to the GoGraph object.
  /*!
    A reference to GO graph to be used.
  */
  std::shared_ptr<const GoGraph> _graph;

  //! The information content map.
  /*!
    An information content map.
  */
  std::shared_ptr<const TermInformationContentMap> _icMap;

  //! A method for calculating the sum of information content of the terms in a set.
  /*!
    This method calculates the sum of information content of the terms in a set.
  */
  [[nodiscard]] double calcICSum(const OntologySetType<std::string> &terms) const {

    double sum{0.0};
    for (auto const &term : terms) {

      sum += _icMap->getValue(term);

    }

    return sum;

  }


};

} // namespace

#endif
