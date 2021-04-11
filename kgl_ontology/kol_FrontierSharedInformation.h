/*=============================================================================
Copyright (c) 2016 Paul W. Bible
Distributed under the Boost Software License, Version 1.0.
==============================================================================*/
#ifndef KGL_FRONTIER_SHARED_INFORMATION
#define KGL_FRONTIER_SHARED_INFORMATION

#include "kol_TermInformationContentMap.h"
#include "SharedInformationInterface.h"
#include "SetUtilities.h"
#include "Accumulators.h"
#include "kol_GoGraph.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kellerberrin::ontology {


/*! \class FrontierSharedInformation
	\brief A class to calculate shared infromation across disjoint common ancestors in linear time.

	This class calculates shared infromation along a semantic frontier between terms.
*/
class FrontierSharedInformation : public SharedInformationInterface {

public:

  //! A constructor
  /*!
    Creates the CoutoGraSMSharedInformation class
  */
  FrontierSharedInformation(const std::shared_ptr<const GoGraph> &goGraph, const std::shared_ptr<const TermInformationContentMap> &icMap)
      : _goGraph(goGraph), _icMap(icMap) {}

  ~FrontierSharedInformation() override = default;


  //! A method for determining the common disjunctive ancestors
  /*!
    This method returns the common disjunctive ancestors for two terms
  */
  OntologySetType<std::string> getCommonDisjointAncestors(const std::string &termC1, const std::string &termC2) const {

    //common disjoint ancestors set
    OntologySetType<std::string> cda;

    //commonDisjointAncestors(c,c) = {c}, by definition
    if (termC1 == termC2) {

      cda.insert(termC1);
      return cda;

    }

    //std::cout << "Linear GraSm " << std::endl;
    OntologySetType<std::string> ancestorsC1 = _goGraph->getSelfAncestorTerms(termC1);
    //std::cout << ancestorsC1.size() << std::endl;
    OntologySetType<std::string> ancestorsC2 = _goGraph->getSelfAncestorTerms(termC2);
    //std::cout << ancestorsC2.size() << std::endl;


    //Couto: Anc = CommonAnc(c1,c2)
    OntologySetType<std::string> commonAncestors = SetUtilities::set_intersection(ancestorsC1, ancestorsC2);
    //std::cout << commonAncestors.size() << std::endl;

    //commonDisjointAncestors(c,c) = {c}, by definition
    if (commonAncestors.size() == 1) {

      return commonAncestors;

    }

    //std::cout << "CA size " << commonAncestors.size() << std::endl;

    //get the boost graph
    const GoGraph::Graph &go_graph = _goGraph->getGraph();

    OntologySetType<std::size_t> edgesC1;
    OntologySetType<std::size_t> edgesC2;

    const GoGraph::EdgeIndexMap &edge_index_map = _goGraph->edgeIndexMap();
    OntologyMapType<std::string, OntologySetType<std::size_t> > termToEdges;

    EdgeSetVisitor c1EdgeVisitor(edgesC1, edge_index_map, termToEdges);
    EdgeSetVisitor c2EdgeVisitor(edgesC2, edge_index_map, termToEdges);

    //get edges for c1
    boost::breadth_first_search(go_graph, _goGraph->getVertexByName(termC1), boost::visitor(c1EdgeVisitor));

    //get edges for c1
    boost::breadth_first_search(go_graph, _goGraph->getVertexByName(termC2), boost::visitor(c2EdgeVisitor));

    //std::cout << "edges 1 " << edgesC1.size() << std::endl;
    //std::cout << "edges 2 " << edgesC2.size() << std::endl;
    //std::cout << "edge map " << termToEdges.size() << std::endl;

    for (auto const &term : commonAncestors) {

      //std::cout << term << std::endl;
      OntologySetType<std::size_t> edges = termToEdges[term];
      //std::cout << edges.size() << std::endl;

      bool isDisj = false;

      //check if the term is attached to a frontier edge
      for (auto const &edgeId : edges) {

        if (edgesC1.find(edgeId) == edgesC1.end() or edgesC2.find(edgeId) == edgesC2.end()) {

          isDisj = true;
          break;

        }

      }

      //if the term is a disjoint ancestor add it to the set
      if (isDisj) {

        cda.insert(term);

      }

    }

    //std::cout << "frontier cda size " << cda.size() << std::endl;
    return cda;

  }

  //! An method for returning the shared information of two terms
  /*!
    This method returns the mean information content of the frontier ancestors
  */
  [[nodiscard]] double sharedInformation(const std::string &termA, const std::string &termB) const override {
    // return 0 for any terms not in the datbase
    if (not _icMap->hasTerm(termA) || not _icMap->hasTerm(termB)) {

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
    if (!_icMap->hasTerm(term)) {

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


  //breadth first search to calculate the visited edges
  class EdgeSetVisitor : public boost::default_bfs_visitor {

  public:
    EdgeSetVisitor(OntologySetType<std::size_t> &inSet,
                   const GoGraph::EdgeIndexMap &inMap,
                   OntologyMapType<std::string, OntologySetType<std::size_t> > &termToEdges) :
        edgeSet(inSet), eMap(inMap), termEdgesMap(termToEdges) {}

    template<typename Edge, typename Graph>
    void examine_edge(const Edge &edge, const Graph &graph) {
      //add the edge to the set of visited edges
      edgeSet.insert(eMap[edge]);

      //get the vertex of the target
      typename Graph::vertex_descriptor v = boost::target(edge, graph);
      std::string term = graph[v].termId;

      //create a set for the term if none exists
      if (termEdgesMap.find(term) == termEdgesMap.end()) {

        termEdgesMap[term] = OntologySetType<std::size_t>();

      }
      //add the edge to the map
      termEdgesMap[term].insert(eMap[edge]);

    }

    OntologySetType<std::size_t> &edgeSet;
    const GoGraph::EdgeIndexMap &eMap;
    OntologyMapType<std::string, OntologySetType<std::size_t> > &termEdgesMap;

  };

  std::shared_ptr<const GoGraph> _goGraph;
  std::shared_ptr<const TermInformationContentMap> _icMap;


};

} // namespace

#endif
