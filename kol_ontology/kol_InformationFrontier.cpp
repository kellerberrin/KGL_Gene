//
// Created by kellerberrin on 16/4/21.
//

#include "kol_InformationFrontier.h"

#include "kol_SetUtilities.h"
#include "kol_Accumulators.h"
#include "kol_GoGraphImpl.h"

#include <utility>
#include <algorithm>

#include <boost/graph/breadth_first_search.hpp>


namespace kol = kellerberrin::ontology;


//breadth first search to calculate the visited edges
class EdgeSetVisitor : public boost::default_bfs_visitor {

public:
  EdgeSetVisitor(kol::OntologySetType<std::size_t> &inSet,
                 const kol::GoGraphImpl::EdgeIndexMap &inMap,
                 kol::OntologyMapType<std::string, kol::OntologySetType<std::size_t> > &termToEdges) :
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

      termEdgesMap[term] = kol::OntologySetType<std::size_t>();

    }
    //add the edge to the map
    termEdgesMap[term].insert(eMap[edge]);

  }

  kol::OntologySetType<std::size_t> &edgeSet;
  const kol::GoGraphImpl::EdgeIndexMap &eMap;
  kol::OntologyMapType<std::string, kol::OntologySetType<std::size_t> > &termEdgesMap;

};


//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/

kol::OntologySetType<std::string> kol::InformationFrontier::getCommonDisjointAncestors(const std::string &termC1,
                                                                                       const std::string &termC2) const {

  //common disjoint ancestors set
  OntologySetType<std::string> cda;

  //commonDisjointAncestors(c,c) = {c}, by definition
  if (termC1 == termC2) {

    cda.insert(termC1);
    return cda;

  }

  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termC1);
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getGoGraphImpl().getSelfAncestorTerms(termC2);

  //Couto: Anc = CommonAnc(c1,c2)
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);

  //commonDisjointAncestors(c,c) = {c}, by definition
  if (commonAncestors.size() == 1) {

    return commonAncestors;

  }


  //get the boost graph
  const GoGraphImpl::Graph &go_graph = graph_ptr_->getGoGraphImpl().getGraph();

  OntologySetType<std::size_t> edgesC1;
  OntologySetType<std::size_t> edgesC2;

  const GoGraphImpl::EdgeIndexMap &edge_index_map = graph_ptr_->getGoGraphImpl().edgeIndexMap();
  OntologyMapType<std::string, OntologySetType<std::size_t> > termToEdges;

  EdgeSetVisitor c1EdgeVisitor(edgesC1, edge_index_map, termToEdges);
  EdgeSetVisitor c2EdgeVisitor(edgesC2, edge_index_map, termToEdges);

  //get edges for c1
  boost::breadth_first_search(go_graph, graph_ptr_->getGoGraphImpl().getVertexByName(termC1), boost::visitor(c1EdgeVisitor));

  //get edges for c1
  boost::breadth_first_search(go_graph, graph_ptr_->getGoGraphImpl().getVertexByName(termC2), boost::visitor(c2EdgeVisitor));


  for (auto const &term : commonAncestors) {

    OntologySetType<std::size_t> edges = termToEdges[term];

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

  return cda;

}

//! An method for returning the shared information of two terms
/*!
  This method returns the mean information content of the frontier ancestors
*/
double kol::InformationFrontier::sharedInformation(const std::string &termA, const std::string &termB) const {
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

