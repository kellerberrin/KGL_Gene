//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_InformationCoutoGraSM.h"

#include "kol_Accumulators.h"
#include "kol_SetUtilities.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>

namespace kol = kellerberrin::ontology;



//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/
kol::OntologySetType<std::string> kol::InformationCoutoGraSM::getCommonDisjointAncestors(const std::string &termC1,
                                                                                         const std::string &termC2) const {

  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getSelfAncestorTerms(termC1);
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getSelfAncestorTerms(termC2);

  //Couto: CommonDisjAnc = {}
  OntologySetType<std::string> cda;

  if (termC1 == termC2) {

    cda.insert(termC1);

    return cda;

  }

  //Couto: Anc = CommonAnc(c1,c2)
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);

  std::vector<std::pair<double, std::string> > orderedCommonAncestors;
  //create a pair to associate a term with its information content
  for (auto const &term : commonAncestors) {

    orderedCommonAncestors.emplace_back(ic_map_ptr_->termInformation(term), term);

  }

  //sort descending
  std::sort(orderedCommonAncestors.begin(), orderedCommonAncestors.end(), std::greater<>());

  //start of main algorithm
  //Couto: for all a in sortDescByIC(Anc) do ...
  for (auto const&[value, termA] : orderedCommonAncestors) {

    //Couto: isDisj=true
    bool isDisj = true;

    //Couto: for all cda in CommonDisjAnc do ...
    for (auto const &termCda : cda) {

      //continue if the terms are the same
      if (termCda == termA) {

        continue;

      }

      //Couto: isDisj = isDisj ^ ( DisjAnc(c1,(a,cda)) or DisjAnc(c2,(a,cda)) )
      isDisj = isDisj && (isDisjoint(termC1, termA, termCda) || isDisjoint(termC2, termA, termCda));

    }

    //Couto: if isDisj then...
    if (isDisj) {

      //Couto: addTo(CommonDisjAnc,a)
      cda.insert(termA);

    }

  }

  return cda;

}


//! A method for determining if for a term c, a pair (a1,a2) is disjoint in c
/*!
  This method returns
*/
bool kol::InformationCoutoGraSM::isDisjoint(const std::string &termC,
                                            const std::string &termA1,
                                            const std::string &termA2) const {

  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(termA1) != graph_ptr_->getTermOntology(termA2) ||
      graph_ptr_->getTermOntology(termC) != graph_ptr_->getTermOntology(termA1) ||
      graph_ptr_->getTermOntology(termC) != graph_ptr_->getTermOntology(termA2)) {

    return false;

  }

  if (ic_map_ptr_->termInformation(termA1) <= ic_map_ptr_->termInformation(termA2)) {

    size_t nPaths = getNumPaths(termA1, termA2);
    size_t nPaths1 = getNumPaths(termA1, termC);
    size_t nPaths2 = getNumPaths(termA2, termC);

    if (nPaths1 >= nPaths * nPaths2) {

      return true;

    } else {

      return false;

    }

  } else {

    return false;

  }

}

//! A method for calculating the number of paths for one term to another.
/*!
  This method returns the number of paths between two terms
*/
size_t kol::InformationCoutoGraSM::getNumPaths(const std::string &termA, const std::string &termB) const {

  if (ic_map_ptr_->termInformation(termA) > ic_map_ptr_->termInformation(termB)) {

    return 0;

  }

  return pathCount(termA, termB);

}

//! An method for returning the shared information of two terms
/*!
  This method returns the mean information content disjoint common ancestors
*/
double kol::InformationCoutoGraSM::sharedInformation(const std::string &termA, const std::string &termB) const {
  // return 0 for any terms not in the datbase
  if (not ic_map_ptr_->validateTerms(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);

  for (auto const &term : cda) {

    meanIC(ic_map_ptr_->termInformation(term));

  }

  return Accumulators::extractMean(meanIC);

}


//! Count paths from B to A
/*!
  Count paths between B and A
*/
size_t kol::InformationCoutoGraSM::pathCount(const std::string &termA, const std::string &termB) const {

  if (ic_map_ptr_->termInformation(termA) > ic_map_ptr_->termInformation(termB)) {

    return 0;

  }

  OntologySetType<std::string> ancestors = graph_ptr_->getAncestorTerms(termB);
  OntologySetType<std::string> finished;
  OntologyMapType<std::string, size_t> pathMap;
  ancestors.insert(termB);
  const GoGraphImpl::Graph &go_graph = graph_ptr_->getGraph();
  GoGraphImpl::GoVertex v = graph_ptr_->getTermRootVertex(termB);
  visitHelper(v, go_graph, ancestors, finished, pathMap);

  return pathMap[termA];

}

//! Recursive helper method that performs the DFS topological sort for path counting
/*!
  A path counting topological sort recursive method.
*/
void kol::InformationCoutoGraSM::visitHelper(const GoGraphImpl::GoVertex &go_vertex,
                                             const GoGraphImpl::Graph &go_graph,
                                             OntologySetType<std::string> &ancestors,
                                             OntologySetType<std::string> &finished,
                                             OntologyMapType<std::string, size_t> &pathMap) const {

  size_t childCount = 0;
  std::string vTerm = go_graph[go_vertex].termId;

  //examine children and recurse
  GoGraphImpl::InEdgeIterator it, end;
  for (boost::tie(it, end) = boost::in_edges(go_vertex, go_graph); it != end; ++it) {

    GoGraphImpl::GoVertex child = boost::source(*it, go_graph);
    std::string childTerm = go_graph[child].termId;
    if (not ancestors.contains(childTerm)) {

      continue;

    }
    //recurse if child is not finished
    if (not finished.contains(childTerm)) {

      visitHelper(child, go_graph, ancestors, finished, pathMap);

    }

    ++childCount;

  }

  //finish vertex
  finished.insert(vTerm);

  if (childCount == 0) {

    pathMap[vTerm] = 1;

  } else {

    pathMap[vTerm] = 0;
    for (boost::tie(it, end) = boost::in_edges(go_vertex, go_graph); it != end; ++it) {
      GoGraphImpl::GoVertex child = boost::source(*it, go_graph);
      std::string childTerm = go_graph[child].termId;

      if (not ancestors.contains(childTerm)) {

        continue;

      }

      pathMap[vTerm] += pathMap[childTerm];

    }

  }

}

//! A private function to create a string key from a pair of terms
/*!
  Creates a string key our of a pair to use in memorizing path counts
*/
std::string kol::InformationCoutoGraSM::keyPair(const std::string &termA, const std::string &termB) const {

  if (termA.compare(termB) > 0) {

    return termB + "_" + termA;

  } else {

    return termB + "_" + termA;

  }

}

//! A private function to test if the key as been seen already
/*!
  A private function to test if the key as been seen already.
*/
bool kol::InformationCoutoGraSM::hasSeenKey(const std::string &key) const {

  if (path_memory_.find(key) != path_memory_.end()) {

    return true;

  } else {

    return false;

  }

}
