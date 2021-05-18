//
// Created by kellerberrin on 19/4/21.
//

#include "kol_OntologyTypes.h"
#include "kol_TermProbabilityMap.h"
#include "kol_Accumulators.h"
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

};//end class, DFS visitor class


// Constructor
void kol::TermProbabilityMap::calcProbabilityMap(const std::shared_ptr<const GoGraph> &graph,
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

  const double const_bp_annotations = getRootCount(GO::getRootTermBP());
  const double const_mf_annotations = getRootCount(GO::getRootTermMF());
  const double const_cc_annotations = getRootCount(GO::getRootTermCC());

  Accumulators::SimpleAccumulator minMaxBP;
  Accumulators::SimpleAccumulator minMaxMF;
  Accumulators::SimpleAccumulator minMaxCC;

  for (auto& [term_id, anno_ont_pair] : probability_map_) {

    auto& [annotations, ontology] = anno_ont_pair;

    if (annotations == 0.0) {

      continue;

    }

    switch (ontology) {

      case GO::Ontology::BIOLOGICAL_PROCESS:
        minMaxBP(annotations);
        annotations = annotations / const_bp_annotations;
        break;

      case GO::Ontology::MOLECULAR_FUNCTION:
        minMaxMF(annotations);
        annotations = annotations / const_mf_annotations;
        break;

      case GO::Ontology::CELLULAR_COMPONENT:
        minMaxCC(annotations);
        annotations = annotations / const_cc_annotations;
        break;

      default:
      case GO::Ontology::ONTO_ERROR:
        ExecEnv::log().error("TermProbabilityMap::calcProbabilityMap; GO term: {} does not have a valid ontology", term_id);
        break;

    }

  }

  //calculate single annotation minimum normalization factors
  bp_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_1anno_ = 1.0 / Accumulators::extractMax(minMaxCC);

  //calculate minimum annotation minimum normalization factors
  bp_normalization_min_min_anno_ = Accumulators::extractMin(minMaxBP) / Accumulators::extractMax(minMaxBP);
  mf_normalization_min_min_anno_ = Accumulators::extractMin(minMaxMF) / Accumulators::extractMax(minMaxMF);
  cc_normalization_min_min_anno_ = Accumulators::extractMin(minMaxCC) / Accumulators::extractMax(minMaxCC);

}//end constructor logic


double kol::TermProbabilityMap::getRootCount(const std::string& root_id) const {

  auto const result = probability_map_.find(root_id);
  if (result == probability_map_.end()) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {} not in probability map", root_id);
    return 1.0;

  }

  auto const& [term_id, anno_ont_pair] = *result;
  auto const& [annotations, ontology] = anno_ont_pair;

  if (annotations <= 0.0) {

    ExecEnv::log().error("TermProbabilityMap::getRootCount; root term: {}, invalid annotation count: {}", root_id, annotations);
    return 1.0;

  }

  return annotations;

}


double kol::TermProbabilityMap::getValue(const std::string &term_id) const {

  auto result = probability_map_.find(term_id);
  if (result == probability_map_.end()) {

    return badIdValue();

  }
  auto const& [map_term_id, prob_ont_pair] = *result;
  auto const& [probability, ontology] = prob_ont_pair;
  return probability;

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


