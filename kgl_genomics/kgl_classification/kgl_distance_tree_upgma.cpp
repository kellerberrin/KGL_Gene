//
// Created by kellerberrin on 16/12/17.
//

#include "kel_exec_env.h"
#include "kgl_distance_tree_upgma.h"

#include <ranges>

namespace kgl = kellerberrin::genome;


void kgl::MatrixGenerator::initializeMatrix(const TreeNodeVector& tree_node_vector) {

  tree_node_vector_ = tree_node_vector;

  calculateMatrix(tree_node_vector_, distance_matrix_);


}

bool kgl::MatrixGenerator::sumMatrix(const TreeNodeVector& tree_node_vector) {

  if (tree_node_vector.size() != tree_node_vector_.size()) {

    ExecEnv::log().warn("Cannot add different node sizes; existing node size: {}, add node size: {}",
                        tree_node_vector_.size(), tree_node_vector.size());
    return false;

  }

  for (auto const& [existing, adding] : std::ranges::views::zip(tree_node_vector_, tree_node_vector)) {

    if (existing->nodeText() != adding->nodeText()) {

      ExecEnv::log().warn("Cannot add dissimilar nodes; existing node: {}, add node: {}",
                          existing->nodeText(), adding->nodeText());
      return false;

    }

  }

  DistanceMatrix add_matrix;
  calculateMatrix(tree_node_vector, add_matrix);

  distance_matrix_.addMatrix(add_matrix);

  return true;

}

void kgl::MatrixGenerator::calculateMatrix(const TreeNodeVector& tree_node_vector, DistanceMatrix& distance_matrix) {

  // Resize.
  distance_matrix.resize(tree_node_vector.size());

  // Populate the matrix.
  for (size_t row = 0; row < tree_node_vector.size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      distance_matrix.setDistance(row, column, distance(tree_node_vector[row], tree_node_vector[column]));

    }

  }

}

kgl::DistanceType_t kgl::MatrixGenerator::distance(const std::shared_ptr<TreeNodeDistance>& row_node,
                                                     const std::shared_ptr<TreeNodeDistance>& column_node) {

  return row_node->distance(column_node);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::DistanceTreeUPGMA::DistanceTreeUPGMA(const MatrixGenerator& tree_matrix) {

  tree_node_vector_ = tree_matrix.treeNodes();
  distance_matrix_.copyMatrix(tree_matrix.distanceMatrix());

}



size_t kgl::DistanceTreeUPGMA::getLeafCount(size_t leaf_idx) const {

  if (leaf_idx >= tree_node_vector_.size()) {

    ExecEnv::log().error("getLeafCount(), bad index: {}, node vector size: {}", leaf_idx, tree_node_vector_.size());
    return 1;

  }

  return tree_node_vector_[leaf_idx]->leafNodeCount();

}


kgl::TreeNodeVector kgl::DistanceTreeUPGMA::calculateTree(bool normalized) {

  if (normalized) {

    distance_matrix_.normalizeDistance();

  }

  UPGMATree();
  return tree_node_vector_;

}


// Reduces the calculateDistance matrix.
// The reduced column is the left most column (column, j index = 0)
void kgl::DistanceTreeUPGMA::reduceDistance(size_t i, size_t j) {

  // Save and resize
  size_t reduce_size = distance_matrix_.size() - 1;

  DistanceMatrix temp_distance(std::move(distance_matrix_));

  distance_matrix_.resize(reduce_size);

  // re-populate merged distances.
  size_t idx_row = 1;
  for(size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      auto i_leaf = static_cast<DistanceType_t>(getLeafCount(i));
      auto j_leaf = static_cast<DistanceType_t>(getLeafCount(j));
      DistanceType_t calc_dist = (i_leaf * temp_distance.getDistance(row, i)) + (j_leaf * temp_distance.getDistance(row, j));
      calc_dist = calc_dist / (i_leaf + j_leaf);
      distance_matrix_.setDistance(idx_row, 0, calc_dist);
      ++idx_row;

    }

  }

  // re-populate other distances.
  if (distance_matrix_.size() <= 2) {

    return;

  }

  idx_row = 2;  // shift down 2 rows.
  bool update = false;
  for (size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      size_t idx_column = 1;  // shift to column = 1
      for (size_t column = 0; column < row; ++column) {

        if (column != i and column != j) {

          distance_matrix_.setDistance(idx_row, idx_column, temp_distance.getDistance(row, column));
          idx_column++;
          update = true;

        }

      }

      if (update) {

        idx_row++;
        update = false;

      }

    }

  }

}


bool kgl::DistanceTreeUPGMA::reduceNode(size_t row, size_t column, DistanceType_t minimum) {

  TreeNodeVector temp_node_vector;
  std::shared_ptr<TreeNodeDistance> column_node = tree_node_vector_[column];
  std::shared_ptr<TreeNodeDistance> row_node = tree_node_vector_[row];

  for (size_t idx = 0; idx < tree_node_vector_.size(); idx++) {

    if (not (idx == column or idx == row)) {

      temp_node_vector.push_back(tree_node_vector_[idx]);

    }

  }

  tree_node_vector_ = std::move(temp_node_vector);

  DistanceType_t node_distance = minimum / 2;
  DistanceType_t row_distance = node_distance - row_node->parentDistance();
  row_node->parentDistance(row_distance);
  DistanceType_t column_distance = node_distance - column_node->parentDistance();
  column_node->parentDistance(column_distance);
  size_t row_leaves = row_node->leafNodeCount();
  size_t column_leaves = column_node->leafNodeCount();
  auto merged_node_ptr = std::make_shared<CladeNode>("Clade Node", row_leaves + column_leaves);
  merged_node_ptr->parentDistance(node_distance);
  merged_node_ptr->addOutNode(row_node);
  row_node->parentNode(merged_node_ptr);
  merged_node_ptr->addOutNode(column_node);
  column_node->parentNode(merged_node_ptr);
  // Insert the merged node at the front of the vector.
  // This matches the pattern of the reduction of the calculateDistance matrix (above).
  tree_node_vector_.insert(tree_node_vector_.begin(), merged_node_ptr);

  return true;

}


void kgl::DistanceTreeUPGMA::UPGMATree() {

  while (tree_node_vector_.size() > 1) {

    auto [min, row, column] = distance_matrix_.minimum();

    reduceDistance(row, column);

    reduceNode(row, column, min);

  } // while reduceNode.

}

