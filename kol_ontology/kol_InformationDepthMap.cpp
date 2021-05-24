//
// Created by kellerberrin on 19/4/21.
//


#include "kol_OntologyTypes.h"
#include "kol_InformationDepthMap.h"

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>


namespace kol = kellerberrin::ontology;


//! Function to return all the keys in the map
/*!
  Returns all valid keys in the map.
*/
std::vector<std::string> kol::InformationDepthMap::getKeys() const {

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
size_t kol::InformationDepthMap::getValue(const std::string &termId) const {

  if (not hasTerm(termId) or termId.empty()) {

    return 0;

  }
  //get index
  auto const&[term, index] = *(name_to_index_.find(termId));
  //return the depth
  return depths_.at(index);

}

//! A method for calculating the least common ancestor
/*!
  This method searches the sets to determine the deepest common ancestor
*/
std::string kol::InformationDepthMap::getLCA(const OntologySetType<std::string> &ancestorsA,
                                             const OntologySetType<std::string> &ancestorsB) const {
  //get the first term as a start
  if (ancestorsA.size() < ancestorsB.size()) {

    return getEfficientLCA(ancestorsA, ancestorsB);

  } else {

    return getEfficientLCA(ancestorsB, ancestorsA);

  }

}

//! A private method to calculate the depth values on object construction
/*!
  This method actually calculates the depth values.
*/

void kol::InformationDepthMap::initializeDepthMap(const GoGraph &graph) {

  //Initialize an annotation list the size of verticies in go, each value is 0
  //depths_ = std::vector<std::size_t>(graph->getNumVertices(),0);
  depths_ = std::vector<size_t>(graph.getNumVertices(), 0);

  //get the (first) root of the ontology.
  //GoGraph::GoVertex root = graph->getRoot();
  //TESTING
  //std::cout << root << std::endl;
  //std::cout << graph->getTermStringIdByIndex(root) << std::endl;

  //get the boost graph from the GoGraph object. Must be done to utilize boost algorithms
  const GoGraph::Graph &go_graph = graph.getGraph();

  //wrap _depth with a vertex map
  const GoGraph::VertexIndexMap &vMap = graph.vertexIndexMap();

  /* // Temparary fix until I can get SWIG to recognize std::size_t
  boost::iterator_property_map< std::vector<std::size_t>::iterator,
                              GoGraph::VertexIndexMap >
                  d_map(depths_.begin(), vMap);
  */
  boost::iterator_property_map<std::vector<size_t>::iterator, GoGraph::VertexIndexMap> d_map(depths_.begin(), vMap);

  //call the boost depth first search using our custom visitor
  // revering the graph is necessary otherwise the root vertex would have no edges.
  //boost::depth_first_search(boost::make_reverse_graph(*go_graph),boost::visitor(vis).root_vertex(root));
  GoGraph::GoVertex bpRoot = graph.getVertexByName(GO::getRootTermBP());
  GoGraph::GoVertex mfRoot = graph.getVertexByName(GO::getRootTermMF());
  GoGraph::GoVertex ccRoot = graph.getVertexByName(GO::getRootTermCC());

  //Start at bproot, record depths
  // must reverse graph due to edge relationship direction
  boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                              bpRoot,
                              boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

  //Start at bproot, record depths
  // must reverse graph due to edge relationship direction
  boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                              mfRoot,
                              boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

  //Start at bproot, record depths
  // must reverse graph due to edge relationship direction
  boost::breadth_first_search(boost::make_reverse_graph(go_graph),
                              ccRoot,
                              boost::visitor(boost::make_bfs_visitor(boost::record_distances(d_map, boost::on_tree_edge()))));

  // Only valid for std::unordered_map.
  // name_to_index_ = OntologyMapType<std::string,std::size_t>(boost::num_vertices(go_graph));

  // Vertex Iterators
  GoGraph::GoVertexIterator vi, vend;
  for (boost::tie(vi, vend) = boost::vertices(go_graph); vi != vend; ++vi) {

    GoGraph::GoVertex v = *vi;
    name_to_index_[graph.getTermStringIdByIndex(vMap[v])] = vMap[v];
    //std::cout << vMap[v] << " " << depths_.at(vMap[v]) << " " << graph->getTermNameByIndex(vMap[v]) << std::endl;
    //std::cin.get();
  }

}


//! A method for calculating the least common ancestor
/*!
  This method searches the sets to determine the deepest common ancestor
*/
std::string kol::InformationDepthMap::getEfficientLCA(const OntologySetType<std::string> &smaller_set, const OntologySetType<std::string> &larger_set) const {
  //get the first term as a start
  std::string lca;
  //max depth
  size_t max{0};

  //loop over shorter list
  for (auto const &currentTerm : smaller_set) {

    if (larger_set.find(currentTerm) != smaller_set.end()) {

      //if new max, update
      if (getValue(currentTerm) > max) {

        lca = currentTerm;
        max = getValue(currentTerm);

      }

    }

  }

  return lca;

}
