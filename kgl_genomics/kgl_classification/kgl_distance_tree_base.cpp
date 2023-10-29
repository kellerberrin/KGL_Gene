//
// Created by kellerberrin on 23/05/19.
//

#include "kgl_distance_tree_base.h"

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VirtualDistanceNode::zeroDistance(std::shared_ptr<const VirtualDistanceNode> node) const {

  if (not node) {

    ExecEnv::log().critical("VirtualDistanceNode::zeroDistance; null pointer passed as argument");

  }

  return dynamic_cast<const void*>(this) == dynamic_cast<const void*>(node.get()) ;

}



// Recursively counts the total number of leaf nodes.
size_t kgl::TreeDistanceNode::leafNodeCount() const {

  size_t leaf_nodes = 0;

  if (not isLeaf()) {

    for (auto node : outNodes()) {

      leaf_nodes += node.second->leafNodeCount();

    }

    return leaf_nodes;

  } else {

    return 1;

  }

}
