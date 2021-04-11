/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_GRASM_SHARED_INFORMATION
#define KGL_GRASM_SHARED_INFORMATION

#include "SharedInformationInterface.h"
#include "kol_TermInformationContentMap.h"
#include "kol_GoGraph.h"
#include "Accumulators.h"
#include "SetUtilities.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kellerberrin::ontology {

/*! \class CoutoGraSMSharedInformation
	\brief A class to calculate shared information across disjoint common ancestors using the exact algorithm as written in the paper.

	This class calculates shared information across disjoint common ancsetors.

    F. M. Couto, M. J. Silva, and P. M. Coutinho, "Measuring semantic similarity
	between Gene Ontology terms," Data & Knowledge Engineering, vol. 61, 
	pp. 137-152, Apr 2007.

	Couto proposed calculating this value as a subsitute for the IC of the MICA in calculating
	 Resnik, Lin, and Jiang-Conrath

*/
class CoutoGraSMSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMGreaterOrEqual class
  */
  CoutoGraSMSharedInformation(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~CoutoGraSMSharedInformation() override = default;

  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  [[nodiscard]] OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const {

    OntologySetType<std::string> ancestorsC1 = _goGraph->getSelfAncestorTerms(termC1);
    //std::cout << ancestorsC1.size() << std::endl;
    OntologySetType<std::string> ancestorsC2 = _goGraph->getSelfAncestorTerms(termC2);
    //std::cout << ancestorsC2.size() << std::endl;

    //Couto: CommonDisjAnc = {}
    OntologySetType<std::string> cda;

    if (termC1 == termC2) {

      cda.insert(termC1);

      return cda;

    }

    //Couto: Anc = CommonAnc(c1,c2)
    OntologySetType<std::string> commonAncestors = SetUtilities::set_intersection(ancestorsC1, ancestorsC2);
    //std::cout << commonAncestors.size() << std::endl;

    std::vector<std::pair<double, std::string> > orderedCommonAncestors;
    //create a pair to associate a term with its information content
    for (auto const &term : commonAncestors) {

      orderedCommonAncestors.emplace_back(_icMap->getValue(term), term);

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


  //! A method for determining if for a term c, a pair (a1,a2) is disjoint in c
  /*!
    This method returns
  */
  [[nodiscard]] bool isDisjoint(const std::string &termC, const std::string &termA1, const std::string &termA2) const {

    //std::cout << "isDisjoint " << termC << " ("  << termA1 << " , " << termA2 << ") "; //<< std::endl;
    //if not from same ontology, return 0;
    if (_goGraph->getTermOntology(termA1) != _goGraph->getTermOntology(termA2) ||
        _goGraph->getTermOntology(termC) != _goGraph->getTermOntology(termA1) ||
        _goGraph->getTermOntology(termC) != _goGraph->getTermOntology(termA2)) {
      return false;
    }

    if (_icMap->getValue(termA1) <= _icMap->getValue(termA2)) {
      //std::cout << "case 1" << std::endl;
      size_t nPaths = getNumPaths(termA1, termA2);
      //std::cout << "nPaths " << termA1 << " to "  << termA2 << " " << nPaths << std::endl << std::endl;
      size_t nPaths1 = getNumPaths(termA1, termC);
      //std::cout << "nPaths " << termA1 << " to "  << termC << " " << nPaths1 << std::endl << std::endl;
      size_t nPaths2 = getNumPaths(termA2, termC);
      //std::cout << "nPaths " << termA2 << " to "  << termC << " " << nPaths2 << std::endl << std::endl;
      if (nPaths1 >= nPaths * nPaths2) {
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


  //! A method for calculating the number of paths for one term to another.
  /*!
    This method returns the number of paths between two terms
  */
  [[nodiscard]] size_t getNumPaths(const std::string &termA, const std::string &termB) const {

    if (_icMap->getValue(termA) > _icMap->getValue(termB)) {

      return 0;

    }

    return pathCount(termA, termB);

  }


  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content disjoint common ancestors
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override {
    // return 0 for any terms not in the datbase
    if (not _icMap->hasTerm(termA) or not _icMap->hasTerm(termB)) {

      return 0.0;

    }
    // return 0 for terms in different ontologies
    if (_goGraph->getTermOntology(termA) != _goGraph->getTermOntology(termB)) {

      return 0.0;

    }

    Accumulators::MeanAccumulator meanIC;
    OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);
    //std::cout << "size " << cda.size() << std::endl;

    for (auto const &term : cda) {
      //std::cout << *iter << std::endl;
      //std::cout << _icMap[*iter] << std::endl;
      meanIC(_icMap->getValue(term));

    }

    return Accumulators::extractMean(meanIC);

  }

  //! An interface method for returning the shared information of a single terms,or information content
  /*!
    This method privdes a mechanism for returing a term's infromation content.
  */
  [[nodiscard]] double sharedInformation(const std::string &term) const override {
    // return 0 for any terms not in the datbase
    if (not _icMap->hasTerm(term)) {

      return 0.0;

    }

    return _icMap->getValue(term);

  }

  //! An interface method for returning the maximum information content for a term
  /*!
    This method provides the absolute max information content within a corpus for normalization purposes.
  */
  [[nodiscard]] double maxInformationContent(const std::string &term) const override {


    //select the correct ontology normalization factor
    GO::Ontology ontology = _goGraph->getTermOntology(term);
    double maxIC;

    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        maxIC = _icMap->getMinBP();
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        maxIC = _icMap->getMinMF();
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        maxIC = _icMap->getMinCC();
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        maxIC = 0.0;
        break;

    }

    if (maxIC <= 0.0) {

      return 0.0;

    }

    return -1.0 * std::log(maxIC);

  }


  //! An interface method for determining if a term can be found
  /*!
    Determines if the term can be found in the current map.
  */
  [[nodiscard]] bool hasTerm(const std::string &term) const override {

    return _icMap->hasTerm(term);

  }

  //! An interface method for determining if the two terms are of like ontologies.
  /*!
    Determine if two terms are of the same ontology.
  */
  [[nodiscard]] bool isSameOntology(const std::string &termA, const std::string &termB) const override {

    return _goGraph->getTermOntology(termA) == _goGraph->getTermOntology(termB);

  }


private:

  //! Count paths from B to A
  /*!
    Count paths between B and A
  */
  [[nodiscard]] size_t pathCount(const std::string &termA, const std::string &termB) const {

    if (_icMap->getValue(termA) > _icMap->getValue(termB)) {

      return 0;

    }

    OntologySetType<std::string> ancestors = _goGraph->getAncestorTerms(termB);
    OntologySetType<std::string> finished;
    OntologyMapType<std::string, size_t> pathMap;
    ancestors.insert(termB);
    const GoGraph::Graph &go_graph = _goGraph->getGraph();
    GoGraph::GoVertex v = _goGraph->getTermRootVertex(termB);
    visitHelper(v, go_graph, ancestors, finished, pathMap);

    return pathMap[termA];

  }

  //! Recursive helper method that performs the DFS topological sort for path counting
  /*!
    A path counting topological sort recursive method.
  */
  void visitHelper(const GoGraph::GoVertex &go_vertex,
                   const GoGraph::Graph &go_graph,
                   OntologySetType<std::string> &ancestors,
                   OntologySetType<std::string> &finished,
                   OntologyMapType<std::string, size_t> &pathMap) const {
    size_t childCount = 0;
    std::string vTerm = go_graph[go_vertex].termId;
    //std::cout << "discover vertex " << vTerm << std::endl;

    //examine children and recurse
    GoGraph::InEdgeIterator it, end;
    for (boost::tie(it, end) = boost::in_edges(go_vertex, go_graph); it != end; ++it) {

      GoGraph::GoVertex child = boost::source(*it, go_graph);
      std::string childTerm = go_graph[child].termId;
      if (!SetUtilities::set_contains(ancestors, childTerm)) {
        continue;
      }
      //recurse if child is not finished
      if (!SetUtilities::set_contains(finished, childTerm)) {

        visitHelper(child, go_graph, ancestors, finished, pathMap);

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
      for (boost::tie(it, end) = boost::in_edges(go_vertex, go_graph); it != end; ++it) {
        GoGraph::GoVertex child = boost::source(*it, go_graph);
        std::string childTerm = go_graph[child].termId;

        if (!SetUtilities::set_contains(ancestors, childTerm)) {
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
  [[nodiscard]] std::string keyPair(const std::string &termA, const std::string &termB) const {

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
  [[nodiscard]] bool hasSeenKey(const std::string &key) const {

    if (_pathMemory.find(key) != _pathMemory.end()) {

      return true;

    } else {

      return false;

    }

  }

  std::shared_ptr<const GoGraph> _goGraph;
  std::shared_ptr<const TermInformationContentMap> _icMap;
  OntologyMapType<std::string, size_t> _pathMemory;

};


} // namespace


#endif
