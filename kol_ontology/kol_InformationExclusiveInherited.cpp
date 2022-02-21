//
// Created by kellerberrin on 19/4/21.
//


#include "kol_InformationExclusiveInherited.h"

#include "kol_SetUtilities.h"
#include "contrib/kol_Accumulators.h"
#include "contrib/kol_GoGraphImpl.h"

#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kol = kellerberrin::ontology;


//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/
kol::OntologySetType<std::string> kol::InformationExclusiveInherited::getCommonDisjointAncestors(const std::string &termC1,
                                                                                                 const std::string &termC2) const {

  //Z&L: EICommonSet <- <empty_set>
  OntologySetType<std::string> cda;

  //EICA(t,t) = {t}
  if (termC1 == termC2) {

    cda.insert(termC1);
    return cda;

  }

  //Z&L: CommonAnSet <- GetCommonAnSet(t1,t2,AncestorSet)
  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termC1);
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termC2);
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);

  //commonDisjointAncestors(c,c) = {c}, by definition
  if (commonAncestors.size() == 1) {

    return commonAncestors;

  }

  //Z&L: UnionAnSet <- GetAnSet(t1,AncestorSet) U GetAnSet(t2,AncestorSet)
  OntologySetType<std::string> unionAncestors = SetUtilities::setUnion(ancestorsC1, ancestorsC2);

  //Z&L: DiffAnSet <- UnionAnSet - CommonAnSet
  OntologySetType<std::string> diffSet = SetUtilities::setDifference(unionAncestors, commonAncestors);

  //get the boost graph
  const GoGraphImpl::Graph &g = graph_ptr_->getGoGraphImpl().getGraph();
  //Z&L: for each a in CommonAnSet do ...
  for (auto const &term : commonAncestors) {

    bool isDisj = false;

    GoGraphImpl::GoVertex v = graph_ptr_->getGoGraphImpl().getVertexByName(term);
    GoGraphImpl::InEdgeIterator ei, end;
    for (boost::tie(ei, end) = boost::in_edges(v, g); ei != end; ++ei) {

      GoGraphImpl::GoVertex child = boost::source(*ei, g);
      size_t index = graph_ptr_->getGoGraphImpl().getVertexIndex(child);
      std::string cTerm = graph_ptr_->getGoGraphImpl().getTermStringIdByIndex(index);

      if (diffSet.find(cTerm) != diffSet.end()) {
        //early exit should improve runtime
        isDisj = true;
        break;

      }

    }


    //if the term is a disjoint ancestor add it to the set
    if (isDisj) {

      cda.insert(term);

    }

  }

  return cda;

}


//! An method for returning the shared information of two terms
/*!
  This method returns the mean information content of the frontier ancestors
*/
double kol::InformationExclusiveInherited::sharedInformation(const std::string &termA, const std::string &termB) const {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->validateTerms(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);

  for (auto const &element : cda) {

    meanIC(ic_map_ptr_->termInformation(element));

  }

  return Accumulators::extractMean(meanIC);

}



