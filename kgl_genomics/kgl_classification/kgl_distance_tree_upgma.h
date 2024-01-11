//
// Created by kellerberrin on 16/12/17.
//

#ifndef KGL_DISTANCE_TREE_UPGMA_H
#define KGL_DISTANCE_TREE_UPGMA_H

#include "kel_utility.h"
#include "kgl_distance_matrix.h"
#include "kgl_distance_tree_node.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate a UPGMA distance tree.
// Initialized with a vector of leaf nodes.
// Generates the tree structure and returns the root of the tree.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceTreeUPGMA {

public:

  DistanceTreeUPGMA() = default;
  ~DistanceTreeUPGMA() = default;

  void addDistanceMap(const TreeNodeVector& tree_node_vector) { tree_node_vector_.clear(); tree_node_vector_ = tree_node_vector; }
  // Returns the root of the calculated tree.
  TreeNodeVector calculateTree();

private:

  TreeNodeVector tree_node_vector_;
  DistanceMatrix distance_matrix_;

  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void reduceDistance(size_t i, size_t j);
  void UPGMATree();

  [[nodiscard]] DistanceType_t distance(const std::shared_ptr<TreeNodeDistance>& row_node, const std::shared_ptr<TreeNodeDistance>& column_node) const;
  [[nodiscard]] size_t getLeafCount(size_t leaf_idx) const;
  void initializeDistance();
  void normalizeDistance();
  void identityZeroDistance();

};


}   // end namespace


#endif //KGL_DISTANCE_TREE_UPGMA_H
