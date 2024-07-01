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
//
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MatrixGenerator {

public:

  MatrixGenerator() = default;
  ~MatrixGenerator() = default;

  bool sumMatrix(const TreeNodeVector& tree_node_vector);
  void initializeMatrix(const TreeNodeVector& tree_node_vector);

  [[nodiscard]] const TreeNodeVector& treeNodes() const { return tree_node_vector_; }
  [[nodiscard]] const DistanceMatrix& distanceMatrix() const { return distance_matrix_; }

private:

  TreeNodeVector tree_node_vector_;
  DistanceMatrix distance_matrix_;

  static void calculateMatrix(const TreeNodeVector& tree_node_vector, DistanceMatrix& distance_matrix);
  [[nodiscard]] static DistanceType_t distance( const std::shared_ptr<TreeNodeDistance>& row_node,
                                                const std::shared_ptr<TreeNodeDistance>& column_node);

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate a UPGMA distance tree.
// Initialized with a vector of leaf nodes.
// Generates the tree structure and returns the root of the tree.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceTreeUPGMA {

public:

  explicit DistanceTreeUPGMA(const MatrixGenerator& tree_matrix);
  ~DistanceTreeUPGMA() =default;

  // Returns the root of the calculated tree. Default scales distances between [0, 1].
  TreeNodeVector calculateTree(bool normalized = true);

private:

  TreeNodeVector tree_node_vector_;
  DistanceMatrix distance_matrix_;

  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void reduceDistance(size_t i, size_t j);
  void UPGMATree();

  [[nodiscard]] size_t getLeafCount(size_t leaf_idx) const;

};


}   // end namespace


#endif //KGL_DISTANCE_TREE_UPGMA_H
