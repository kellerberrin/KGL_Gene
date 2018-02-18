//
// Created by kellerberrin on 16/12/17.
//


#include <boost/numeric/ublas/triangular.hpp>
// #include <boost/numeric/ublas/io.hpp>

#include "kgl_exec_env.h"
#include "kgl_statistics_upgma.h"


namespace kgl = kellerberrin::genome;
namespace bnu = boost::numeric::ublas;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation of a strict lower triangular distance matrix using the Boost library.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using DistanceMatrixImplType = bnu::triangular_matrix<kgl::DistanceType_t, bnu::strict_lower>;

class kgl::DistanceMatrix::BoostDistanceMatrix {

public:

  explicit BoostDistanceMatrix(size_t matrix_size) : lower_triangular_(matrix_size, matrix_size) {}
  explicit BoostDistanceMatrix(const BoostDistanceMatrix&) = default;
  ~BoostDistanceMatrix() = default;

  size_t size() const { return lower_triangular_.size1(); }
  void resize(size_t new_size) { lower_triangular_.resize(new_size, new_size, false); }
  kgl::DistanceType_t getDistance(size_t i, size_t j) const;
  void setDistance(size_t i, size_t j, kgl::DistanceType_t distance);

private:

  DistanceMatrixImplType lower_triangular_;

};


kgl::DistanceType_t kgl::DistanceMatrix::BoostDistanceMatrix::getDistance(size_t i, size_t j) const {

  if (i >= size() or i >= size()) {

    kgl::ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
    return 0;

  }

  if (i == j) {

    return 0;

  }

  if (j > i) {

    std::swap(i, j);

  }

  return lower_triangular_(i,j);

}

void kgl::DistanceMatrix::BoostDistanceMatrix::setDistance(size_t i, size_t j, kgl::DistanceType_t distance) {

  if (i >= size() or i >= size()) {

    kgl::ExecEnv::log().error("Index too large i: {}, j:{} for distance matrix size: {}", i, j, size());
    return;

  }

  if (i == j) {

    return;

  }

  if (j > i) {

    std::swap(i, j);

  }

  lower_triangular_(i,j) = distance;

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Classification Nodes.
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public class of the UPGMA distance matrix.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

kgl::DistanceMatrix::DistanceMatrix(size_t matrix_size) : diagonal_impl_ptr_(std::make_unique<kgl::DistanceMatrix::BoostDistanceMatrix>(matrix_size)) {}

kgl::DistanceMatrix::DistanceMatrix(const DistanceMatrix& copy) : diagonal_impl_ptr_(std::make_unique<kgl::DistanceMatrix::BoostDistanceMatrix>(*copy.diagonal_impl_ptr_)) {}

kgl::DistanceMatrix::~DistanceMatrix() {}  // DO NOT DELETE or USE DEFAULT. Required because of incomplete PIMPL type.


size_t kgl::DistanceMatrix::size() const {

  return diagonal_impl_ptr_->size();

}


void kgl::DistanceMatrix::resize(size_t new_size) {

  diagonal_impl_ptr_->resize(new_size);

}


kgl::DistanceType_t kgl::DistanceMatrix::getDistance(size_t i, size_t j) const {

  return diagonal_impl_ptr_->getDistance( i, j);

}

void kgl::DistanceMatrix::setDistance(size_t i, size_t j, DistanceType_t value) {

  diagonal_impl_ptr_->setDistance( i, j, value);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance Implementation functions.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::DistanceMatrix::minimum(size_t& i, size_t& j) const {

  bool first_pass = true;
  kgl::DistanceType_t min_distance = 0;
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        min_distance = getDistance( row, column);
        i = row;
        j = column;
        first_pass = false;

      } else {

        kgl::DistanceType_t check_min = getDistance( row, column);
        if (check_min < min_distance) {

          min_distance = check_min;
          i = row;
          j = column;

        }

      }

    }

  }

  return min_distance;

}




kgl::DistanceType_t kgl::DistanceMatrix::maximum(size_t& i, size_t& j) const {

  bool first_pass = true;
  kgl::DistanceType_t max_distance = 0;
  for (size_t row = 0; row < size(); ++row) {

    for (size_t column = 0; column < row; ++column) {

      if (first_pass) {

        max_distance = getDistance( row, column);
        i = row;
        j = column;
        first_pass = false;

      } else {

        kgl::DistanceType_t check_max = getDistance( row, column);
        if (check_max > max_distance) {

          max_distance = check_max;
          i = row;
          j = column;

        }

      }

    }

  }

  return max_distance;

}


// Reduces the distance matrix.
// The reduced column is the left most column (column, j index = 0)
void kgl::DistanceMatrix::reduce(size_t i, size_t j) {

  // Save and resize
  DistanceMatrix temp_distance(*this);

  size_t reduce_size = size() - 1;

  resize(reduce_size);

  // re-populate merged distances.
  size_t idx_row = 1;
  for(size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      kgl::DistanceType_t i_leaf = getLeafCount(i);
      kgl::DistanceType_t j_leaf = getLeafCount(j);
      kgl::DistanceType_t calc_dist = (i_leaf * temp_distance.getDistance(row, i)) + (j_leaf * temp_distance.getDistance(row, j));
      calc_dist = calc_dist / (i_leaf + j_leaf);
      setDistance(idx_row, 0, calc_dist);
      ++idx_row;

    }

  }

  // re-populate other distances.
  if (size() <= 2) {

    return;

  }

  idx_row = 2;  // shift down 2 rows.
  bool update = false;
  for (size_t row = 0; row < temp_distance.size(); ++row) {

    if (row != i and row != j) {

      size_t idx_column = 1;  // shift to column = 1
      for (size_t column = 0; column < row; ++column) {

        if (column != i and column != j) {

          setDistance(idx_row, idx_column, temp_distance.getDistance(row, column));
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t kgl::UPGMAMatrix::getLeafCount(size_t leaf_idx) const {

  if (leaf_idx >= node_vector_ptr_->size()) {

    ExecEnv::log().error("getLeafCount(), bad index: {}, node vector size: {}", leaf_idx, node_vector_ptr_->size());
    return 1;
  }

  return node_vector_ptr_->at(leaf_idx)->leafNodeCount();

}


kgl::DistanceType_t kgl::UPGMAMatrix::distance(std::shared_ptr<PhyloNode> row_node,
                                     std::shared_ptr<PhyloNode> column_node) const {

  return row_node->leaf()->distance(column_node->leaf());

}


kgl::DistanceType_t kgl::UPGMAMatrix::zeroDistance(std::shared_ptr<PhyloNode> row_node,
                                                   std::shared_ptr<PhyloNode> column_node) const {

  return row_node->leaf()->zeroDistance(column_node->leaf());

}


void kgl::UPGMAMatrix::initializeDistance() {

  for (size_t row = 0; row < node_vector_ptr_->size(); ++row) {

    for (size_t column = 0; column < row; column++) {

      setDistance(row, column, distance(node_vector_ptr_->at(row), node_vector_ptr_->at(column)));

    }

  }

  normalizeDistance();

}


void kgl::UPGMAMatrix::normalizeDistance() {

  rescaleDistance();
  identityZeroDistance();

}



void kgl::UPGMAMatrix::identityZeroDistance() {

  for (size_t row = 0; row < node_vector_ptr_->size(); ++row) {

    for (size_t column = 0; column < row; column++) {

      if (zeroDistance(node_vector_ptr_->at(row), node_vector_ptr_->at(column)) == 0) {

          setDistance(row, column, 0.0);

      }

    }

  }

}


void kgl::UPGMAMatrix::rescaleDistance() {

  size_t row;
  size_t column;
  DistanceType_t min = minimum(row, column);
  DistanceType_t max = maximum(row, column);
  DistanceType_t range = max - min;

  if (range == 0.0) {

    ExecEnv::log().error("UPGMAMatrix::normalizeDistance() distance range for all nodes is zero");
    return;

  }

  for (size_t row = 0; row < node_vector_ptr_->size(); ++row) {

    for (size_t column = 0; column < row; column++) {

      DistanceType_t raw_distance = getDistance(row, column);
      DistanceType_t adj_distance = (raw_distance - min) / range;
      setDistance(row, column, adj_distance);

    }

  }

}



void kgl::UPGMAMatrix::calculateReduce() {

  while (node_vector_ptr_->size() > 1) {

    size_t row;
    size_t column;
    DistanceType_t min = minimum(row, column);

    reduce(row, column);

    reduceNode(row, column, min);

  } // while reduceNode.

}


bool kgl::UPGMAMatrix::reduceNode(size_t row, size_t column, DistanceType_t minimum) {

  std::shared_ptr<PhyloNodeVector> temp_node_vector_ptr(std::make_shared<PhyloNodeVector>());
  std::shared_ptr<PhyloNode> column_node = nullptr;
  std::shared_ptr<PhyloNode> row_node = nullptr;

  for (size_t idx = 0; idx < node_vector_ptr_->size(); idx++) {

    if (not (idx == column or idx == row)) {

      temp_node_vector_ptr->push_back(node_vector_ptr_->at(idx));

    } else if (idx == column) {

      column_node = node_vector_ptr_->at(idx);

    } else if (idx == row) {

      row_node = node_vector_ptr_->at(idx);

    }

  }

  if (not column_node or not row_node) {

    ExecEnv::log().error("Null pointer found, error calculating UPGMA, row: {}, column: {}, nodes: {}",
                         row, column, node_vector_ptr_->size());
    return false;

  }

  node_vector_ptr_ = temp_node_vector_ptr;

  DistanceType_t node_distance = minimum / 2;
  DistanceType_t row_distance = node_distance - row_node->distance();
  row_node->distance(row_distance);
  DistanceType_t column_distance = node_distance - column_node->distance();
  column_node->distance(column_distance);
  std::shared_ptr<PhyloNode> merged_node(std::make_shared<PhyloNode>(row_node->leaf()));
  merged_node->distance(node_distance);
  merged_node->addOutNode(row_node);
  merged_node->addOutNode(column_node);
  // Insert the merged node at the front of the vector.
  // This matches the pattern of the reduction of the distance matrix (above).
  node_vector_ptr_->insert(node_vector_ptr_->begin(), merged_node);

  return true;

}


bool kgl::UPGMAMatrix::writeNewick(const std::string& file_name) const {

  std::ofstream newick_file;

  // Open input file.

  newick_file.open(file_name);

  if (not newick_file.good()) {

    ExecEnv::log().error("I/O error; could not open Newick file: {}", file_name);
    return false;

  }

  for (auto node : *node_vector_ptr_) {

    writeNode(node, newick_file);

  }

  newick_file << ";";

  newick_file.close();

  return true;

}


void kgl::UPGMAMatrix::writeNode(std::shared_ptr<PhyloNode> node, std::ofstream& newick_file) const {

  if (node->getMap().size() > 0) {

    newick_file << "(";

    bool first_pass = true;
    for (auto child_node : node->getMap()) {

      if (first_pass) {

        first_pass = false;

      } else {

        newick_file << ",";

      }

      writeNode(child_node.second, newick_file);

    }

    newick_file << ")";

  } else {

    node->leaf()->writeNode(newick_file);

  }

  newick_file << ":";
  newick_file << node->distance();

}

