//
// kgl_distance_tree_node.cpp — tree classification nodes.
//
// Created by kellerberrin on 28/12/23.
//

#include "kgl_distance_tree_node.h"

#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization level namespace


DistanceType_t TreeNodeDistance::distance(const std::shared_ptr<const TreeNodeDistance>&) const {

  ExecEnv::log().warn("Probable logic error; default node distance called.");
  return 0;

}


// Recursively counts the total number of leaf nodes.
size_t TreeNodeDistance::leafNodeCount() const {

  if (isLeaf()) {

    return 1;

  }

  size_t leaf_nodes = 0;
  for (auto const& [distance, out_node_ptr] : outNodes()) {

    leaf_nodes += out_node_ptr->leafNodeCount();

  }

  return leaf_nodes;

}


}   // end namespace