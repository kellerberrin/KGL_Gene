//
// Created by kellerberrin on 16/12/17.
//

#include "kel_exec_env.h"
#include "kgl_distance_tree_upgma.h"


namespace kgl = kellerberrin::genome;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t kgl::DistanceTreeUPGMA::getLeafCount(size_t leaf_idx) const {

  if (leaf_idx >= node_vector_ptr_->size()) {

    ExecEnv::log().error("getLeafCount(), bad index: {}, node vector size: {}", leaf_idx, node_vector_ptr_->size());
    return 1;

  }

  return node_vector_ptr_->at(leaf_idx)->leafNodeCount();

}


kgl::DistanceType_t kgl::DistanceTreeUPGMA::distance(std::shared_ptr<TreeDistanceNode> row_node,
                                                     std::shared_ptr<TreeDistanceNode> column_node) const {

  return row_node->node()->distance(column_node->node());

}


void kgl::DistanceTreeUPGMA::initializeDistance() {

  for (size_t row = 0; row < node_vector_ptr_->size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      distance_matrix_.setDistance(row, column, distance(node_vector_ptr_->at(row), node_vector_ptr_->at(column)));

    }

  }

  normalizeDistance();

}


void kgl::DistanceTreeUPGMA::normalizeDistance() {

  distance_matrix_.normalizeDistance();
  identityZeroDistance();

}

// Enforce the distance metric; d(a, a) = 0.
void kgl::DistanceTreeUPGMA::identityZeroDistance() {

  for (size_t row = 0; row < node_vector_ptr_->size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      auto row_hash = std::hash<std::shared_ptr<const VirtualDistanceNode>>()(node_vector_ptr_->at(row)->node());
      auto column_hash = std::hash<std::shared_ptr<const VirtualDistanceNode>>()(node_vector_ptr_->at(column)->node());

      if (row_hash == column_hash) {

        distance_matrix_.setDistance(row, column, 0.0);

      }

    }

  }

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

  std::shared_ptr<DistanceNodeVector> temp_node_vector_ptr(std::make_shared<DistanceNodeVector>());
  std::shared_ptr<TreeDistanceNode> column_node = node_vector_ptr_->at(column);
  std::shared_ptr<TreeDistanceNode> row_node = node_vector_ptr_->at(row);

  for (size_t idx = 0; idx < node_vector_ptr_->size(); idx++) {

    if (not (idx == column or idx == row)) {

      temp_node_vector_ptr->push_back(node_vector_ptr_->at(idx));

    }

  }

  node_vector_ptr_ = temp_node_vector_ptr;

  DistanceType_t node_distance = minimum / 2;
  DistanceType_t row_distance = node_distance - row_node->distance();
  row_node->distance(row_distance);
  DistanceType_t column_distance = node_distance - column_node->distance();
  column_node->distance(column_distance);
  std::shared_ptr<TreeDistanceNode> merged_node(std::make_shared<TreeDistanceNode>(row_node->node()));
  merged_node->distance(node_distance);
  merged_node->addOutNode(row_node);
  merged_node->addOutNode(column_node);
  // Insert the merged node at the front of the vector.
  // This matches the pattern of the reduction of the calculateDistance matrix (above).
  node_vector_ptr_->insert(node_vector_ptr_->begin(), merged_node);

  return true;

}


void kgl::DistanceTreeUPGMA::UPGMATree() {

  while (node_vector_ptr_->size() > 1) {

    auto [min, row, column] = distance_matrix_.minimum();

    reduceDistance(row, column);

    reduceNode(row, column, min);

  } // while reduceNode.

}


bool kgl::DistanceTreeUPGMA::writeNewick(const std::string& file_name) const {

  std::ofstream newick_file;

  newick_file.open(file_name);
  if (not newick_file.good()) {

    ExecEnv::log().error("I/O error; could not open Newick file: {}", file_name);
    return false;

  }

  for (const auto& node : *node_vector_ptr_) {

    writeNode(node, newick_file);

  }

  newick_file << ";";

  return true;

}


void kgl::DistanceTreeUPGMA::writeNode(const std::shared_ptr<TreeDistanceNode>& node, std::ofstream& newick_file) const {

  if (not node->outNodes().empty()) {

    newick_file << "(";

    bool first_pass = true;
    for (const auto& [distance, child_node] : node->outNodes()) {

      if (first_pass) {

        first_pass = false;

      } else {

        newick_file << ",";

      }

      writeNode(child_node, newick_file);

    }

    newick_file << ")";

  } else {

    node->node()->writeNode(newick_file);

  }

  newick_file << ":";
  newick_file << node->distance();

}

