//
// Created by kellerberrin on 16/4/21.
//

#include "kol_FrontierSharedInformation.h"

namespace kol = kellerberrin::ontology;


//breadth first search to calculate the visited edges
class EdgeSetVisitor : public boost::default_bfs_visitor {

public:
  EdgeSetVisitor(kol::OntologySetType<std::size_t> &inSet,
                 const kol::GoGraph::EdgeIndexMap &inMap,
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
  const kol::GoGraph::EdgeIndexMap &eMap;
  kol::OntologyMapType<std::string, kol::OntologySetType<std::size_t> > &termEdgesMap;

};


//! A method for determining the common disjunctive ancestors
/*!
  This method returns the common disjunctive ancestors for two terms
*/

kol::OntologySetType<std::string> kol::FrontierSharedInformation::getCommonDisjointAncestors(const std::string &termC1,
                                                                                             const std::string &termC2) const {

  //common disjoint ancestors set
  OntologySetType<std::string> cda;

  //commonDisjointAncestors(c,c) = {c}, by definition
  if (termC1 == termC2) {

    cda.insert(termC1);
    return cda;

  }

  //std::cout << "Linear GraSm " << std::endl;
  OntologySetType<std::string> ancestorsC1 = graph_ptr_->getSelfAncestorTerms(termC1);
  //std::cout << ancestorsC1.size() << std::endl;
  OntologySetType<std::string> ancestorsC2 = graph_ptr_->getSelfAncestorTerms(termC2);
  //std::cout << ancestorsC2.size() << std::endl;


  //Couto: Anc = CommonAnc(c1,c2)
  OntologySetType<std::string> commonAncestors = SetUtilities::setIntersection(ancestorsC1, ancestorsC2);
  //std::cout << commonAncestors.size() << std::endl;

  //commonDisjointAncestors(c,c) = {c}, by definition
  if (commonAncestors.size() == 1) {

    return commonAncestors;

  }

  //std::cout << "CA size " << commonAncestors.size() << std::endl;

  //get the boost graph
  const GoGraph::Graph &go_graph = graph_ptr_->getGraph();

  OntologySetType<std::size_t> edgesC1;
  OntologySetType<std::size_t> edgesC2;

  const GoGraph::EdgeIndexMap &edge_index_map = graph_ptr_->edgeIndexMap();
  OntologyMapType<std::string, OntologySetType<std::size_t> > termToEdges;

  EdgeSetVisitor c1EdgeVisitor(edgesC1, edge_index_map, termToEdges);
  EdgeSetVisitor c2EdgeVisitor(edgesC2, edge_index_map, termToEdges);

  //get edges for c1
  boost::breadth_first_search(go_graph, graph_ptr_->getVertexByName(termC1), boost::visitor(c1EdgeVisitor));

  //get edges for c1
  boost::breadth_first_search(go_graph, graph_ptr_->getVertexByName(termC2), boost::visitor(c2EdgeVisitor));

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
double kol::FrontierSharedInformation::sharedInformation(const std::string &termA, const std::string &termB) const {
// return 0 for any terms not in the datbase
  if (not ic_map_ptr_->hasTerm(termA) || not ic_map_ptr_->hasTerm(termB)) {

    return 0.0;

  }
// return 0 for terms in different ontologies
  if (graph_ptr_->getTermOntology(termA) != graph_ptr_->getTermOntology(termB)) {

    return 0.0;

  }

  Accumulators::MeanAccumulator meanIC;
  OntologySetType<std::string> cda = getCommonDisjointAncestors(termA, termB);
//std::cout << "size " << cda.size() << std::endl;

  for (auto const &term : cda) {
//std::cout << ic_map_ptr_[*iter] << std::endl;
    meanIC(ic_map_ptr_->getValue(term));

  }

  return Accumulators::extractMean(meanIC);

}

//! An interface method for returning the shared information of a single terms,or information content
/*!
  This method privdes a mechanism for returing a term's infromation content.
*/

double kol::FrontierSharedInformation::sharedInformation(const std::string &term) const {
// return 0 for any terms not in the datbase
  if (!ic_map_ptr_->hasTerm(term)) {

    return 0.0;

  }

  return ic_map_ptr_->getValue(term);

}

//! An interface method for returning the maximum information content for a term
/*!
  This method provides the absolute max information content within a corpus for normalization purposes.
*/
double kol::FrontierSharedInformation::maxInformationContent(const std::string &term) const {


//select the correct ontology normalization factor
  GO::Ontology ontology = graph_ptr_->getTermOntology(term);
  double maxIC;

  switch (ontology) {

    case GO::Ontology::BIOLOGICAL_PROCESS:
      maxIC = ic_map_ptr_->getMinBP();
      break;

    case GO::Ontology::MOLECULAR_FUNCTION:
      maxIC = ic_map_ptr_->getMinMF();
      break;

    case GO::Ontology::CELLULAR_COMPONENT:
      maxIC = ic_map_ptr_->getMinCC();
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
