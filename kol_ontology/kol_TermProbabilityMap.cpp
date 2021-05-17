//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_TermProbabilityMap.h"
#include "kel_exec_env.h"

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/reverse_graph.hpp>

namespace kol = kellerberrin::ontology;
namespace kel = kellerberrin;

//! Depth first search boost visitor
/*!
  This defines a class used by TermProbabilityMap to calculate the cumulative annotations
   of a term based on the true path rule. Basically this class is passed to a boost
   Depth first search algorithm to add the number of a child's annotaitons to the parent.
*/

class CumulativeAnnotationsVisitor : public boost::default_dfs_visitor {

public:
  //! A parameterized constructor passing parameters to the boost default_dfs_visitor
  CumulativeAnnotationsVisitor(const std::shared_ptr<const kol::GoGraph> &graph_ptr,
                               const std::shared_ptr<const kol::AnnotationData> &anno_data_ptr,
                               std::vector<std::size_t>& annotations,
                               kol::OntologyMapType<std::string, std::size_t>& name_to_index)
      : graph_ptr_(graph_ptr),
        anno_data_ptr_(anno_data_ptr),
        annotations_(annotations),
        name_index_map_(name_to_index) {}

  //! The extended method of the default_dfs_visitor, finish_vertex
  /*!
    This method is called during a depth first search traversal of a graph and called
      when the visitor finished or leaves a node. This is the last time the algorithm touches
      the node.

      Here we calculate the number of cummulative annotations for a term (node) by adding the current
        term's annotations to the annotations of the term's children.
  */
  template<typename Vertex, typename Graph>
  void finish_vertex(Vertex vertex, const Graph &graph) {

    //Get the vertx index.
    size_t vertex_index = get(boost::vertex_index, graph)[vertex];

    //Get the term_id string.
    std::string term_id = graph_ptr_->getTermStringIdByIndex(vertex_index);

    //Use the AnnotationData object to get the actual number of annotations.
    size_t term_annotations = anno_data_ptr_->getNumAnnotationsForGoTerm(term_id);

    //initialized this map for the TermProbablilityMap, can initialize it for free here
    auto result = name_index_map_.try_emplace(term_id, vertex_index);
    if (not result.second) {

      kel::ExecEnv::log().error("CumulativeAnnotationsVisitor::finish_vertex; unable to add duplicate go term: {} to name index map", term_id);

    }

    //Iterate over the children of the term to add the child annotations
    auto [edge_iter, edge_end] = boost::out_edges(vertex, graph);
    for(;edge_iter != edge_end; ++edge_iter) {

      Vertex child_vertex = boost::target(*edge_iter, graph);

      size_t child_vertex_index = get(boost::vertex_index, graph)[child_vertex];

      std::string child_term_id = graph_ptr_->getTermStringIdByIndex(child_vertex_index);

      term_annotations += annotations_[child_vertex_index];

    } //end for, each child

    annotations_[vertex_index] = term_annotations;

  }//end method, finish_vertex

  //! The go graph object
  std::shared_ptr<const kol::GoGraph> graph_ptr_;

  //! An AnnotationData object for accessing annotations
  std::shared_ptr<const kol::AnnotationData> anno_data_ptr_;

  //! The annotation list to hold the and query the cummulative annotations
  std::vector<std::size_t>& annotations_;

  //! A map from name To index, initialized in this visitor.
  kol::OntologyMapType<std::string, std::size_t>& name_index_map_;

};//end class, DFS visitor class


// Constructor
void kol::TermProbabilityMap::calcProbabilityIndex( const std::shared_ptr<const GoGraph> &graph,
                                                    const std::shared_ptr<const AnnotationData> &annotation_data) {

  // Set the default minimum probablity policy for normalization
  is_single_anno_min_ = false;

  // Initialize an annotation list the size of verticies in go, each value is 0
  probabilities_ = std::vector<double>(graph->getNumVertices(), 0);

  // Get the (first) root of the ontology.
  GoGraph::GoVertex root = graph->getRoot();

  // A variable for the cumulative annotations of the graph.
  std::vector<std::size_t> annotationCounts(graph->getNumVertices(), 0);

  // Create the visitor object.
  CumulativeAnnotationsVisitor visitor(graph, annotation_data, annotationCounts, name_to_index_);

  // Get the boost graph from the GoGraph object. Must be done to utilize boost algorithms.
  const GoGraph::Graph &go_graph = graph->getGraph();

  // Reversing the graph is necessary otherwise the root vertex would have no edges.
  auto reversed_graph = boost::make_reverse_graph(go_graph);
  // Call the boost depth first search using our custom visitor.
  boost::depth_first_search(reversed_graph, boost::visitor(visitor).root_vertex(root));

  // Traverse the vertices to find the roots for each ontology
  Accumulators::SimpleAccumulator minMaxBP;
  Accumulators::SimpleAccumulator minMaxMF;
  Accumulators::SimpleAccumulator minMaxCC;

  // Iterate through all vertices in the graph
  auto [vert_iter, vert_end] = boost::vertices(go_graph);
  for(;vert_iter != vert_end; ++vert_iter) {

    GoGraph::GoVertex vertex = *vert_iter;
    size_t vertex_index = graph->getVertexIndex(vertex);

    auto counts = static_cast<double>(annotationCounts[vertex_index]);
    if (counts == 0.0) {

      continue;

    }

    //switch on ontology, find the max for each set
    switch (graph->getTermOntologyByIndex(vertex_index)) {

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
  auto [vertex_iter, vertex_end] = boost::vertices(go_graph);
  for(;vertex_iter != vertex_end; ++vertex_iter) {

    GoGraph::GoVertex vertex = *vertex_iter;
    std::size_t vertex_index = graph->getVertexIndex(vertex);

    switch (graph->getTermOntologyByIndex(vertex_index)) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        probabilities_[vertex_index] = static_cast<double>(annotationCounts[vertex_index]) / Accumulators::extractMax(minMaxBP);
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        probabilities_[vertex_index] = static_cast<double>(annotationCounts[vertex_index]) / Accumulators::extractMax(minMaxMF);
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        probabilities_[vertex_index] = static_cast<double>(annotationCounts[vertex_index]) / Accumulators::extractMax(minMaxCC);
        break;

      case GO::Ontology::ONTO_ERROR:
        break;

    }

  }

}//end constructor logic


class AnnotationMapVisitor : public boost::default_dfs_visitor {

public:
  //! A parameterized constructor passing parameters to the boost default_dfs_visitor
  AnnotationMapVisitor( const std::shared_ptr<const kol::GoGraph> &graph_ptr,
                        const std::shared_ptr<const kol::AnnotationData> &anno_data_ptr,
                        kol::TermProbOntMap& probability_map)
      : graph_ptr_(graph_ptr),
        anno_data_ptr_(anno_data_ptr),
        probability_map_(probability_map) {}

  //! The extended method of the default_dfs_visitor, finish_vertex
  /*!
    This method is called during a depth first search traversal of a graph and called
      when the visitor finished or leaves a node. This is the last time the algorithm touches
      the node.

      Here we calculate the number of cummulative annotations for a term (node) by adding the current
        term's annotations to the annotations of the term's children.
  */
  template<typename Vertex, typename Graph>
  void finish_vertex(Vertex vertex, const Graph &graph) {

    //Get the vertx index.
    size_t vertex_index = get(boost::vertex_index, graph)[vertex];

    //Get the term_id string.
    std::string term_id = graph_ptr_->getTermStringIdByIndex(vertex_index);

    auto term_ontology = graph_ptr_->getTermOntologyByIndex(vertex_index);
    //Use the AnnotationData object to get the actual number of annotations.
    size_t term_annotations = anno_data_ptr_->getNumAnnotationsForGoTerm(term_id);
    auto float_annotations = static_cast<double>(term_annotations);

    //Iterate over the children of the term to add the child annotations
    auto [edge_iter, edge_end] = boost::out_edges(vertex, graph);
    for(;edge_iter != edge_end; ++edge_iter) {

      Vertex child_vertex = boost::target(*edge_iter, graph);

      size_t child_vertex_index = get(boost::vertex_index, graph)[child_vertex];

      std::string child_term_id = graph_ptr_->getTermStringIdByIndex(child_vertex_index);

      auto child_result = counted_terms_.find(child_term_id);
      if (child_result != counted_terms_.end()) {

        continue;

      }

      auto find_result = probability_map_.find(child_term_id);
      if (find_result == probability_map_.end()) {

        kel::ExecEnv::log().error("CumulativeAnnotationsVisitor::finish_vertex; probability map cannot find term: {}", child_term_id);

      } else {

        auto const& [term_id, anno_ont_pair] = *find_result;
        auto const& [child_annotations, child_ontology] = anno_ont_pair;
        if (child_ontology == term_ontology) {

          float_annotations += child_annotations;

        }

      }
      counted_terms_.insert(child_term_id);

    } //end for, each child

    auto value_pair = std::pair<double, kol::GO::Ontology>{float_annotations, term_ontology};
    auto map_result = probability_map_.try_emplace(term_id, value_pair);
    if (not map_result.second) {

      kel::ExecEnv::log().error("CumulativeAnnotationsVisitor::finish_vertex; unable to add duplicate go term: {} to probability map", term_id);

    }

  }//end method, finish_vertex

  //! The go graph object
  std::shared_ptr<const kol::GoGraph> graph_ptr_;

  //! An AnnotationData object for accessing annotations
  std::shared_ptr<const kol::AnnotationData> anno_data_ptr_;

  //! A map from name To probabilities.
  kol::TermProbOntMap& probability_map_;

  //! Set of counted terms.
  kol::OntologySetType<std::string> counted_terms_;

};//end class, DFS visitor class



// Constructor
void kol::TermProbabilityMap::calcProbabilityMap( const std::shared_ptr<const GoGraph> &graph,
                                                  const std::shared_ptr<const AnnotationData> &annotation_data) {


  GoGraph::GoVertex root = graph->getRoot();

  // Create the visitor object.
  AnnotationMapVisitor map_visitor(graph, annotation_data, probability_map_);

  // Get the boost graph from the GoGraph object. Must be done to utilize boost algorithms.
  const GoGraph::Graph &go_graph = graph->getGraph();

  // Reversing the graph is necessary otherwise the root vertex would have no edges.
  auto reversed_graph = boost::make_reverse_graph(go_graph);
  // Call the boost depth first search using our custom visitor.
  boost::depth_first_search(reversed_graph, boost::visitor(map_visitor).root_vertex(root));

  auto const bp_result = probability_map_.find(GO::getRootTermBP());
  if (bp_result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::calcProbabilityMap; BP root term: {} not in probability map", GO::getRootTermBP());
    return;

  }
  auto const& [bp_term_id, bp_anno_ont_pair] = *bp_result;
  auto const& [bp_annotations, bp_ontology] = bp_anno_ont_pair;
  const double const_bp_annotations{bp_annotations};

  ExecEnv::log().info("TermProbabilityMap::calcProbabilityMap; BP root term: {}, annotation count: {}", bp_term_id, bp_annotations);

  auto const mf_result = probability_map_.find(GO::getRootTermMF());
  if (mf_result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::calcProbabilityMap; MF root term: {} not in probability map", GO::getRootTermMF());
    return;

  }
  auto const& [mf_term_id, mf_anno_ont_pair] = *mf_result;
  auto const& [mf_annotations, mf_ontology] = mf_anno_ont_pair;
  const double const_mf_annotations{mf_annotations};

  ExecEnv::log().info("TermProbabilityMap::calcProbabilityMap; MF root term: {}, annotation count: {}", mf_term_id, mf_annotations);

  auto const cc_result = probability_map_.find(GO::getRootTermCC());
  if (cc_result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::calcProbabilityMap; CC root term: {} not in probability map", GO::getRootTermCC());
    return;

  }
  auto const& [cc_term_id, cc_anno_ont_pair] = *cc_result;
  auto const& [cc_annotations, cc_ontology] = cc_anno_ont_pair;
  const double const_cc_annotations{cc_annotations};

  ExecEnv::log().info("TermProbabilityMap::calcProbabilityMap; CC root term: {}, annotation count: {}", cc_term_id, cc_annotations);

  for (auto& [term_id, anno_ont_pair] : probability_map_) {

    auto& [annotations, ontology] = anno_ont_pair;

    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        annotations = annotations / const_bp_annotations;
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        annotations = annotations / const_mf_annotations;
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        annotations = annotations / const_cc_annotations;
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        ExecEnv::log().error("TermProbabilityMap::calcProbabilityMap; GO term: {} does not have a valid ontology", term_id);
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
