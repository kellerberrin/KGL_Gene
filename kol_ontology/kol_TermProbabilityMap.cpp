//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_TermProbabilityMap.h"

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

namespace kol = kellerberrin::ontology;

//! Depth first search boost visitor
/*!
  This defines a class used by TermProabilityMap to calculate the cumulative annotations
   of a term based on the true path rule. Basically this class is passed to a boost
   Depth first search algorithm to add the number of a child's annotaitons to the parent.
*/

class dfs_cumulative_annotations_visitor : public boost::default_dfs_visitor {

public:
  //! A parameterized constructor passing parameters to the boost default_dfs_visitor
  dfs_cumulative_annotations_visitor(const std::shared_ptr<const kol::GoGraph> &inGraph,
                                     const std::shared_ptr<const kol::AnnotationData> &inData,
                                     std::vector<std::size_t> *annotations,
                                     kol::OntologyMapType<std::string, std::size_t> *nameToIndex)
      : goGraph(inGraph), annoData(inData), annoList(annotations), nameIndexMap(nameToIndex) {}

  //! The extended method of the default_dfs_visitor, finish_vertex
  /*!
    This method is called during a depth first search traversal of a graph and called
      when the visitor finished or leaves a node. This is the last time the algorithm touches
      the node.

      Here we calculate the number of cummulative annotations for a term (node) by adding the current
        term's annotations to the annotations of the term's children.
  */
  template<typename Vertex, typename Graph>
  void finish_vertex(Vertex u, const Graph &g) {

    //get the vertx index
    std::size_t uIndex = get(boost::vertex_index, g)[u];

    //get the termId string
    std::string termId = goGraph->getTermStringIdByIndex(uIndex);

    //use the AnnotationData object to get the actual number of annotations
    std::size_t currentTermAnnos = annoData->getNumAnnotationsForGoTerm(termId);

    //std::cout << termId << std::endl;

    //initialized this map for the TermProbablilityMap, can initialize it for free here
    nameIndexMap->insert(std::make_pair(termId, uIndex));

    //Iterate over the children of the term to add the child annotations
    typename boost::graph_traits<Graph>::out_edge_iterator ei, e_end;
    for (tie(ei, e_end) = boost::out_edges(u, g); ei != e_end; ++ei) {

      Vertex v = boost::target(*ei, g);
      std::size_t vIndex = get(boost::vertex_index, g)[v];

      std::string childTermId = goGraph->getTermStringIdByIndex(vIndex);

      //std::cout << "child " << childTermId << std::endl;

      currentTermAnnos += annoList->at(vIndex);

    }//end for, each child

    annoList->at(uIndex) = currentTermAnnos;

    //std::cout << uIndex << std::endl;
    //std::cout << termId << std::endl;
    //std::cout << annoList->at(uIndex) << std::endl;
    //std::cin.get();

  }//end method, finish_vertex

  //! The go graph object
  std::shared_ptr<const kol::GoGraph> goGraph;

  //! An AnnotationData object for accessing annotations
  std::shared_ptr<const kol::AnnotationData> annoData;

  //! The annotaiton list to hold the and query the cummulative annotations
  std::vector<std::size_t> *annoList;

  //! A map from name To index, initialized in this visitor.
  kol::OntologyMapType<std::string, std::size_t> *nameIndexMap;

};//end class, DFS visitor class


// Constructor
kol::TermProbabilityMap::TermProbabilityMap( const std::shared_ptr<const GoGraph> &graph,
                                             const std::shared_ptr<const AnnotationData> &annoData) {

  //set the default minimum probablity policy for normalization
  is_single_anno_min_ = false;

  //Initialize an annotation list the size of verticies in go, each value is 0
  probabilities_ = std::vector<double>(graph->getNumVertices(), 0);

  //get the (first) root of the ontology.
  GoGraph::GoVertex root = graph->getRoot();

  //TESTING
  //std::cout << root << std::endl;
  //std::cout << graph->getTermStringIdByIndex(root) << std::endl;

  //a variable for the cumulative annotations of the graph
  std::vector<std::size_t> annotationCounts(graph->getNumVertices(), 0);

  //create the visitor object
  dfs_cumulative_annotations_visitor vis(graph, annoData, &annotationCounts, &name_to_index_);

  //get the boost graph from the GoGraph object. Must be done to utilize boost algorithms
  const GoGraph::Graph &go_graph = graph->getGraph();

  //call the boost depth first search using our custom visitor
  // reversing the graph is necessary otherwise the root vertex would have no edges.
  boost::depth_first_search(boost::make_reverse_graph(go_graph), boost::visitor(vis).root_vertex(root));

  // Traverse the vertices to find the roots for each ontology
  Accumulators::SimpleAccumulator minMaxBP;
  Accumulators::SimpleAccumulator minMaxMF;
  Accumulators::SimpleAccumulator minMaxCC;

  // Vertex Iterators
  GoGraph::GoVertexIterator vi, vend;
  for (boost::tie(vi, vend) = boost::vertices(go_graph); vi != vend; ++vi) {

    GoGraph::GoVertex v = *vi;
    std::size_t index = graph->getVertexIndex(v);

    auto counts = static_cast<double>(annotationCounts.at(index));
    if (counts == 0.0) {

      continue;

    }

    //switch on ontology, find the max for each set
    switch (graph->getTermOntologyByIndex(index)) {
      case GO::Ontology::BIOLOGICAL_PROCESS:
        minMaxBP(counts);
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        minMaxMF(counts);
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        minMaxCC(counts);
        break;

      case GO::Ontology::ONTO_ERROR:
        break;
    }

  }//end for, vertex iterator

  //calculate single annotation minimum normalization factors
  bp_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxCC);

  //calculate minimum annotation minimum normalization factors
  bp_normalization_min_min_anno_ = Accumulators::extractMin(minMaxBP) / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_min_anno_ = Accumulators::extractMin(minMaxMF) / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_min_anno_ = Accumulators::extractMin(minMaxCC) / Accumulators::extractMax(minMaxCC);

  //Traverse the cummulative annotation vector, convert to probability, add to probability vector
  for (boost::tie(vi, vend) = boost::vertices(go_graph); vi != vend; ++vi) {

    GoGraph::GoVertex v = *vi;
    std::size_t index = graph->getVertexIndex(v);

    switch (graph->getTermOntologyByIndex(index)) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        probabilities_.at(index) = annotationCounts.at(index) / Accumulators::extractMax(minMaxBP);
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        probabilities_.at(index) = annotationCounts.at(index) / Accumulators::extractMax(minMaxMF);
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        probabilities_.at(index) = annotationCounts.at(index) / Accumulators::extractMax(minMaxCC);
        break;

      case GO::Ontology::ONTO_ERROR:
        break;

    }

  }

}//end constructor logic


//! Function to return all the keys in the map
/*!
  Returns all valid keys in the map.
*/
std::vector<std::string> kol::TermProbabilityMap::getKeys() const {

  std::vector<std::string> keys;
  for (auto const&[key, value] : name_to_index_) {

    keys.push_back(key);

  }

  return keys;

}

//! Mapping function to return the value mapped by key
/*!
Get the value mapped by the given key. A specified function for the [] operator
*/
double kol::TermProbabilityMap::getValue(const std::string &termId) const {

  if (not hasTerm(termId)) {

    return badIdValue();

  }

  auto const&[term, index] = *(name_to_index_.find(termId));

  return probabilities_.at(index);

}

//! Public method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/

double kol::TermProbabilityMap::getMICAinfo( const OntologySetType<std::string> &ancestorsA,
                                             const OntologySetType<std::string> &ancestorsB) const {

  if (ancestorsA.empty() or ancestorsB.empty()) {

    return 0.0;

  }

  // Choose the smaller and larger set for maximum efficiency
  if (ancestorsA.size() < ancestorsB.size()) {

    return getEfficientMICA(ancestorsA, ancestorsB);

  } else {

    return getEfficientMICA(ancestorsB, ancestorsA);

  }

}

//! Private method for calculating the most informative common ancestor value
/*!
  This method searches the sets to determine the most informative ancestor.
*/
double kol::TermProbabilityMap::getEfficientMICA( const OntologySetType<std::string> &smaller_set,
                                                  const OntologySetType<std::string> &larger_set) const {

  double max{0.0};
  //loop over shorter list
  for (auto const &term : smaller_set) {

    if (larger_set.find(term) != larger_set.end()) {

      double term_value = getValue(term);
      if (term_value > max) {

        max = term_value;

      }

    }

  }

  return max;

}
