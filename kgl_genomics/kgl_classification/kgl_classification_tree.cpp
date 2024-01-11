//
// Created by kellerberrin on 23/05/19.
//

#include "kgl_classification_tree.h"

namespace kgl = kellerberrin::genome;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ClassificationTree::writeNewick(const std::string& file_name, size_t max_depth) const {

  std::ofstream newick_file;

  newick_file.open(file_name);
  if (not newick_file.good()) {

    ExecEnv::log().error("I/O error; could not open Newick file: {}", file_name);
    return false;

  }

  if (tree_root_vector_.size() > 1) {


    ExecEnv::log().warn("Writing newick file: {}; expected root nodes = 1,  actual root nodes : {}",
                        file_name, tree_root_vector_.size());

  }

  size_t current_depth{0};
  for (const auto& node : tree_root_vector_) {

    writeNode(node, newick_file, max_depth, current_depth);

  }

  newick_file << ";";

  return true;

}

// Recursive.
void kgl::ClassificationTree::writeNode(const std::shared_ptr<TreeNodeDistance>& node,
                                        std::ofstream& newick_file,
                                        size_t max_depth,
                                        size_t current_depth) const {


  ++current_depth;
  if (current_depth >= max_depth) {

    if (node->isBranch()) {

      newick_file << "Clade_Depth_" << current_depth << "_Leaves_" << node->leafNodeCount();

    } else {

      newick_file << node->nodeText();

    }

  } else if (not node->outNodes().empty()) {

    newick_file << "(";

    bool first_pass = true;
    for (const auto& [distance, child_node] : node->outNodes()) {

      if (first_pass) {

        first_pass = false;

      } else {

        newick_file << ",";

      }

      // Recursion
      writeNode(child_node, newick_file, max_depth, current_depth);

    }

    newick_file << ")";

  } else {

    newick_file << node->nodeText();

  }

  newick_file << ":";
  newick_file << node->parentDistance();

}

