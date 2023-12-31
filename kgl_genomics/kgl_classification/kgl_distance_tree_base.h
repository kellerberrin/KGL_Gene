//
// Created by kellerberrin on 23/05/19.
//

#ifndef KGL_DISTANCE_TREE_H
#define KGL_DISTANCE_TREE_H


#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kel_exec_env.h"
#include "kgl_distance_tree_node.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract Phylogenetic Tree
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceTreeBase {

public:

  explicit DistanceTreeBase(TreeNodeVector tree_root_vector) : tree_root_vector_(std::move(tree_root_vector)) {}
  ~DistanceTreeBase() = default;

  [[nodiscard]] bool writeNewick(const std::string &file_name, size_t max_depth = MAX_TREE_DEPTH_) const;

protected:

  // Top level tree vector.
  // If size() == 1 then this is a rooted tree.
  TreeNodeVector tree_root_vector_;

  constexpr static const size_t MINIMUM_TREE_NODES_{4};
  constexpr static const size_t MAX_TREE_DEPTH_{std::numeric_limits<size_t>::max()};

  void writeNode(const std::shared_ptr<TreeNodeDistance>& node,
                 std::ofstream& newick_file,
                 size_t max_depth,
                 size_t current_depth) const;

};




}   // end namespace


#endif //KGL_DISTANCE_TREE_H
