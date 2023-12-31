//
// Created by kellerberrin on 28/12/23.
//

#include "kgl_distance_tree_node.h"

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::TreeNodeDistance::distance(const std::shared_ptr<const TreeNodeDistance>&) const {

  ExecEnv::log().warn("Probable logic error; default node distance called.");
  return 0;

}


// Recursively counts the total number of leaf nodes.
size_t kgl::TreeNodeDistance::leafNodeCount() const {

  size_t leaf_nodes = 0;

  if (not isLeaf()) {

    for (auto const& [distance, out_node_ptr] : outNodes()) {

      leaf_nodes += out_node_ptr->leafNodeCount();

    }

    return leaf_nodes;

  } else {

    return 1;

  }

}


