//
// Created by kellerberrin on 23/05/19.
//

#include "kgl_distance_tree_base.h"

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Recursively counts the total number of leaf nodes.
size_t kgl::TreeDistanceNode::leafNodeCount() const {

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
