//
// kgl_distance_tree_upgma.cpp — UPGMA distance tree.
//
// Created by kellerberrin on 16/12/17.
//

#include "kel_exec_env.h"
#include "kgl_distance_tree_upgma.h"

#include <ranges>


namespace kellerberrin::genome {   //  organization level namespace


void MatrixGenerator::initializeMatrix(const TreeNodeVector& tree_node_vector) {

  tree_node_vector_ = tree_node_vector;
  calculateMatrix(tree_node_vector_, distance_matrix_);

}


bool MatrixGenerator::sumMatrix(const TreeNodeVector& tree_node_vector) {

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


// Populates the matrix with pairwise distances between all nodes.
void MatrixGenerator::calculateMatrix(const TreeNodeVector& tree_node_vector, DistanceMatrix& distance_matrix) {

  // Resize.
  distance_matrix.resize(tree_node_vector.size());

  // Populate the matrix.
  for (size_t row = 0; row < tree_node_vector.size(); ++row) {
    for (size_t column = 0; column < row; ++column) {

      distance_matrix.setDistance(row, column, tree_node_vector[row]->distance(tree_node_vector[column]));

    }
  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


DistanceTreeUPGMA::DistanceTreeUPGMA(const MatrixGenerator& tree_matrix) {

  tree_node_vector_ = tree_matrix.treeNodes();
  distance_matrix_.copyMatrix(tree_matrix.distanceMatrix());

}


size_t DistanceTreeUPGMA::getLeafCount(size_t leaf_idx) const {

  if (leaf_idx >= tree_node_vector_.size()) {

    ExecEnv::log().error("getLeafCount(), bad index: {}, node vector size: {}", leaf_idx, tree_node_vector_.size());
    return 1;

  }

  return tree_node_vector_[leaf_idx]->leafNodeCount();

}


TreeNodeVector DistanceTreeUPGMA::calculateTree(bool normalized) {

  if (normalized) {

    distance_matrix_.normalizeDistance();

  }

  UPGMATree();
  return tree_node_vector_;

}


// Reduces the distance matrix after merging nodes i and j.
// The merged node occupies the left most column (column index = 0) and the first row.
// The reference implementation performed this reduction with an incremental index
// dance (idx_row starting at 2 with an update flag); this rewrite uses an explicit
// old->new index remap which produces the identical reduced matrix.
void DistanceTreeUPGMA::reduceDistance(size_t i, size_t j) {

  // The merged (row, column) pair must be removed from the matrix.
  const size_t reduce_size = distance_matrix_.size() - 1;

  // The UPGMA weighted average of the merged pair's distances to every surviving node.
  const auto i_leaf = static_cast<DistanceType_t>(getLeafCount(i));
  const auto j_leaf = static_cast<DistanceType_t>(getLeafCount(j));

  // Build the surviving old indices in increasing order.
  // Old survivors map to new indices 1, 2, ...; new index 0 is the merged node.
  std::vector<size_t> new_to_old;
  new_to_old.reserve(reduce_size);
  for (size_t idx = 0; idx < distance_matrix_.size(); ++idx) {

    if (idx != i and idx != j) {

      new_to_old.push_back(idx);

    }

  }

  // Save and resize.
  DistanceMatrix temp_distance(std::move(distance_matrix_));
  distance_matrix_.resize(reduce_size);

  // Re-populate merged distances: new row m (m >= 1), column 0.
  for (size_t new_row = 1; new_row < reduce_size; ++new_row) {

    size_t row = new_to_old[new_row - 1];
    DistanceType_t calc_dist = (i_leaf * temp_distance.getDistance(row, i)) + (j_leaf * temp_distance.getDistance(row, j));
    calc_dist = calc_dist / (i_leaf + j_leaf);
    distance_matrix_.setDistance(new_row, 0, calc_dist);

  }

  // Re-populate the other distances: new [a+1][b+1] = old [survivor_a][survivor_b], a > b.
  if (reduce_size > 2) {

    for (size_t a = 1; a < new_to_old.size(); ++a) {
      for (size_t b = 0; b < a; ++b) {

        distance_matrix_.setDistance(a + 1, b + 1, temp_distance.getDistance(new_to_old[a], new_to_old[b]));

      }
    }

  }

}


// Merges the pair of nodes (i, j) into a new clade node; the vector of tree nodes is updated
// to match the reduced matrix pattern: merged node first, then the survivors in order.
bool DistanceTreeUPGMA::reduceNode(size_t row, size_t column, DistanceType_t minimum) {

  TreeNodeVector temp_node_vector;
  std::shared_ptr<TreeNodeDistance> column_node = tree_node_vector_[column];
  std::shared_ptr<TreeNodeDistance> row_node = tree_node_vector_[row];

  for (size_t idx = 0; idx < tree_node_vector_.size(); idx++) {

    if (not (idx == column or idx == row)) {

      temp_node_vector.push_back(tree_node_vector_[idx]);

    }

  }

  tree_node_vector_ = std::move(temp_node_vector);

  // UPGMA parent-distance bookkeeping. The parent distance of a node doubles as its
  // "height" (distance from the leaves) until the node is attached to a parent.
  DistanceType_t node_distance = minimum / 2;
  DistanceType_t row_distance = node_distance - row_node->parentDistance();
  row_node->parentDistance(row_distance);
  DistanceType_t column_distance = node_distance - column_node->parentDistance();
  column_node->parentDistance(column_distance);
  if (row_distance < 0.0 or column_distance < 0.0) {

    ExecEnv::log().warn("UPGMA negative branch length for node: {} ({}) or node: {} ({})",
                        row_node->nodeText(), row_distance, column_node->nodeText(), column_distance);

  }
  size_t row_leaves = row_node->leafNodeCount();
  size_t column_leaves = column_node->leafNodeCount();
  auto merged_node_ptr = std::make_shared<CladeNode>("Clade Node", row_leaves + column_leaves);
  merged_node_ptr->parentDistance(node_distance);
  merged_node_ptr->addOutNode(row_node);
  row_node->parentNode(merged_node_ptr);
  merged_node_ptr->addOutNode(column_node);
  column_node->parentNode(merged_node_ptr);
  // Insert the merged node at the front of the vector.
  // This matches the pattern of the reduction of the distance matrix (above).
  tree_node_vector_.insert(tree_node_vector_.begin(), merged_node_ptr);

  return true;

}


void DistanceTreeUPGMA::UPGMATree() {

  while (tree_node_vector_.size() > 1) {

    auto [min, row, column] = distance_matrix_.minimum();

    reduceDistance(row, column);

    reduceNode(row, column, min);

  } // while reduceNode.

}


}   // end namespace