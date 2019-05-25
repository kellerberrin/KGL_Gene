//
// Created by kellerberrin on 23/05/19.
//

#include "kgl_phylogenetic_tree.h"

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Recursively counts the total number of leaf nodes.
size_t kgl::PhyloNode::leafNodeCount() const {

  size_t leaf_nodes = 0;

  if (not isLeaf()) {

    for (auto node : getMap()) {

      leaf_nodes += node.second->leafNodeCount();

    }

    return leaf_nodes;

  } else {

    return 1;

  }

}
