//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_ExclusivelyInheritedSharedInformation.h"

namespace kol = kellerberrin::ontology;


//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/
kol::OntologySetType<std::string> kol::ExclusivelyInheritedSharedInformation::getCommonDisjointAncestors( const std::string &termC1,
                                                                                                          const std::string &termC2) const {

  //Z&L: EICommonSet <- <empty_set>
  OntologySetType<std::string> cda;

  //EICA(t,t) = {t}
  if (termC1 == termC2) {

    cda.insert(termC1);
    return cda;

  }

  //Z&L: CommonAnSet <- GetCommonAnSet(t1,t2,AncestorSet)
  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getSelfAncestorTerms(termC1);
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getSelfAncestorTerms(termC2);
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
  const GoGraph::Graph &g = graph_ptr_->getGraph();
  //Z&L: for each a in CommonAnSet do ...
  for (auto const &term : commonAncestors) {

    bool isDisj = false;

    GoGraph::GoVertex v = graph_ptr_->getVertexByName(term);
    GoGraph::InEdgeIterator ei, end;
    for (boost::tie(ei, end) = boost::in_edges(v, g); ei != end; ++ei) {

      GoGraph::GoVertex child = boost::source(*ei, g);
      size_t index = graph_ptr_->getVertexIndex(child);
      std::string cTerm = graph_ptr_->getTermStringIdByIndex(index);

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
double kol::ExclusivelyInheritedSharedInformation::sharedInformation(const std::string &termA, const std::string &termB) const {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->validateTerms(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);

  for (auto const &element : cda) {

    meanIC(ic_map_ptr_->getValue(element));

  }

  return Accumulators::extractMean(meanIC);

}


//! An interface method for returning the shared information of a single terms,or information content
/*!
  This method privdes a mechanism for returing a term's infromation content.
*/

double kol::ExclusivelyInheritedSharedInformation::sharedInformation(const std::string &term) const  {

  return ic_map_ptr_->getValue(term);

}


//! An interface method for returning the maximum information content for a term
/*!
  This method provides the absolute max information content within a corpus for normalization purposes.
*/
double kol::ExclusivelyInheritedSharedInformation::maxInformationContent(const std::string &term) const  {

  return ic_map_ptr_->getMaxInformation(term);

}

