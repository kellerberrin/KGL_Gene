//
// Created by kellerberrin on 19/4/21.
//


#include "kol_InformationCoutoGraSMAdjusted.h"

#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"

#include <boost/graph/breadth_first_search.hpp>


#include <utility>
#include <algorithm>



namespace kol = kellerberrin::ontology;

//! Calculate disjunctive ancestors.
/*!
  A method for determining common disjunctive ancestors for two terms
*/

kol::OntologySetType<std::string> kol::InformationCoutoGraSMAdjusted::getCommonDisjointAncestors(const std::string &termC1,
                                                                                                 const std::string &termC2) const {

  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getAncestorTerms(termC1);
  ancestorsC1.insert(termC1);
  //std::cout << ancestorsC1.size() << std::endl;
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getAncestorTerms(termC2);
  ancestorsC2.insert(termC2);
  //std::cout << ancestorsC2.size() << std::endl;

  //Couto: CommonDisjAnc = {}
  OntologySetType<std::string> cda;

  if (termC1.compare(termC2) == 0) {
    cda.insert(termC1);
    return cda;
  }

  //Couto: Anc = CommonAnc(c1,c2)
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);
  //std::cout << commonAncestors.size() << std::endl;

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

    //std::cout << "testing " << termA << std::endl;

    //Couto: for all cda in CommonDisjAnc do ...
    for (auto const &termCda : cda) {

      //std::cout << "VS " << termCda << std::endl;

      //continue if the terms are the same
      if (termCda == termA) {

        continue;

      }

      //Couto: isDisj = isDisj ^ ( DisjAnc(c1,(a,cda)) or DisjAnc(c2,(a,cda)) )
      isDisj = isDisj && (isDisjoint(termC1, termA, termCda) || isDisjoint(termC2, termA, termCda));

    }

    //Couto: if isDisj then...
    if (isDisj) {
      //std::cout << myPair.second << " is cda " << std::endl;
      //Couto: addTo(CommonDisjAnc,a)
      cda.insert(termA);

    }

  }

  return cda;

}


//! Determine if a terms are disjoint in a concept.
/*!
  A method for determining if, for a term c, a pair (a1,a2) is disjoint in c
*/
bool kol::InformationCoutoGraSMAdjusted::isDisjoint(const std::string &termC,
                                                    const std::string &termA1,
                                                    const std::string &termA2) const {

  //std::cout << "isDisjoint " << termC << " ("  << termA1 << " , " << termA2 << ") "; //<< std::endl;
  //if not from same ontology, return 0;
  if (graph_ptr_->getTermOntology(termA1) != graph_ptr_->getTermOntology(termA2) ||
      graph_ptr_->getTermOntology(termC) != graph_ptr_->getTermOntology(termA1) ||
      graph_ptr_->getTermOntology(termC) != graph_ptr_->getTermOntology(termA2)) {

    return false;

  }

  if (ic_map_ptr_->termInformation(termA1) <= ic_map_ptr_->termInformation(termA2)) {
    //std::cout << "case 1" << std::endl;
    size_t nPaths = getNumPaths(termA1, termA2);
    //std::cout << "nPaths " << termA1 << " to "  << termA2 << " " << nPaths << std::endl << std::endl;
    size_t nPaths1 = getNumPaths(termA1, termC);
    //std::cout << "nPaths " << termA1 << " to "  << termC << " " << nPaths1 << std::endl << std::endl;
    size_t nPaths2 = getNumPaths(termA2, termC);
    //std::cout << "nPaths " << termA2 << " to "  << termC << " " << nPaths2 << std::endl << std::endl;
    if (nPaths1 > nPaths * nPaths2) {
      //std::cout << "true" << std::endl;
      return true;

    } else {
      //std::cout << "false" << std::endl;
      return false;

    }
    //return nPaths1 > nPaths*nPaths2;
  } else {

    return false;

  }

}

//! Calculate the number of paths between two concept terms.
/*!
  A method for calculating the number of paths from one term to another.
*/
size_t kol::InformationCoutoGraSMAdjusted::getNumPaths(const std::string &termA,
                                                       const std::string &termB) const {

  if (ic_map_ptr_->termInformation(termA) > ic_map_ptr_->termInformation(termB)) {

    return 0;

  }

  return pathCount(termA, termB);

}


//! Shared infromation between two conecepts.
/*!
  A method for calculating the shared infromation between two concepts.
*/
double kol::InformationCoutoGraSMAdjusted::sharedInformation(const std::string &termA,
                                                             const std::string &termB) const {
  // return 0 for any terms not in the datbase
  if (!ic_map_ptr_->validateTerms(termA, termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);
  //std::cout << "size " << cda.size() << std::endl;

  for (auto const &term : cda) {
    //std::cout << ic_map_ptr_[*iter] << std::endl;
    meanIC(ic_map_ptr_->termInformation(term));

  }

  return Accumulators::extractMean(meanIC);

}


//! Count paths from B to A
/*!
  Count paths between B and A
*/
std::size_t kol::InformationCoutoGraSMAdjusted::pathCount(const std::string &termA, const std::string &termB) const {

  if (ic_map_ptr_->termInformation(termA) > ic_map_ptr_->termInformation(termB)) {

    return 0;

  }

  OntologySetType<std::string> ancestors = graph_ptr_->getAncestorTerms(termB);
  ancestors.insert(termB);

  const GoGraph::Graph &graph = graph_ptr_->getGraph();
  GoGraph::GoVertex root_vertex = graph_ptr_->getTermRootVertex(termB);

  OntologySetType<std::string> finished;
  OntologyMapType<std::string, size_t> pathMap;
  visitHelper(root_vertex, graph, ancestors, finished, pathMap);

  return pathMap[termA];

}

//! Recursive helper method that performs the DFS topological sort for path counting
/*!
  A path counting topological sort recursive method.
*/
void kol::InformationCoutoGraSMAdjusted::visitHelper(const GoGraph::GoVertex &v,
                                                     const GoGraph::Graph &graph,
                                                     OntologySetType<std::string> &ancestors,
                                                     OntologySetType<std::string> &finished,
                                                     OntologyMapType<std::string, size_t> &pathMap) const {
  size_t childCount = 0;
  std::string vTerm = graph[v].termId;
  //std::cout << "discover vertex " << vTerm << std::endl;

  //examine children and recurse
  GoGraph::InEdgeIterator it, end;
  for (boost::tie(it, end) = boost::in_edges(v, graph); it != end; ++it) {

    GoGraph::GoVertex child = boost::source(*it, graph);
    std::string childTerm = graph[child].termId;

    if (not ancestors.contains(childTerm)) {

      continue;

    }
    //recurse if child is not finished
    if (not finished.contains(childTerm)) {

      visitHelper(child, graph, ancestors, finished, pathMap);

    }
    ++childCount;
  }

  //finish vertex
  finished.insert(vTerm);
  //std::cout << "finish vertex " << vTerm << ", childred " << childCount << std::endl;
  if (childCount == 0) {

    pathMap[vTerm] = 1;

  } else {

    pathMap[vTerm] = 0;
    for (boost::tie(it, end) = boost::in_edges(v, graph); it != end; ++it) {

      GoGraph::GoVertex child = boost::source(*it, graph);
      std::string childTerm = graph[child].termId;
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
std::string kol::InformationCoutoGraSMAdjusted::keyPair(const std::string &termA, const std::string &termB) const {

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
bool kol::InformationCoutoGraSMAdjusted::hasSeenKey(const std::string &key) const {

  if (path_memory_.find(key) != path_memory_.end()) {

    return true;

  } else {

    return false;

  }

}
