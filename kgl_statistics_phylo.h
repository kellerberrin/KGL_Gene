//
// Created by kellerberrin on 6/12/17.
//

#ifndef KGL_STATISTICS_PHYLO_H
#define KGL_STATISTICS_PHYLO_H


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
//using DistanceType_t = uint64_t;
using DistanceType_t = double;

template<class T> class PhyloNode;  // fwd.
template <class T> using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<PhyloNode<T>>>;

template<class T> class PhyloNode {

public:

  explicit PhyloNode(std::shared_ptr<T> leaf, DistanceType_t distance) : leaf_(leaf), distance_(distance) {}
  ~PhyloNode() = default;

  DistanceType_t distance() const { return distance_; }
  void addNode(std::shared_ptr<PhyloNode<T>> node) {
    out_nodes_.insert(std::pair<DistanceType_t , std::shared_ptr<PhyloNode<T>>>(node->distance(), node));
  }
  std::shared_ptr<T> leaf() const { return leaf_; }
  const OutNodes<T>& getMap() const { return out_nodes_; }

private:

  std::shared_ptr<T> leaf_;
  DistanceType_t distance_;
  OutNodes<T> out_nodes_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern to isolate Boost functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceMatrix {

public:

  explicit DistanceMatrix(size_t matrix_size);
  ~DistanceMatrix();  // De not use the default destructor, see PIMPL fwd decl below.

  DistanceType_t getDistance(size_t i, size_t j) const;
  void setDistance(size_t i, size_t j, DistanceType_t distance);
  size_t size() const;

  DistanceType_t reduceMinimum(size_t& i, size_t& j) const;

private:

  class DistanceMatrixImpl;       // Forward declaration of the boost strict diagonal implementation class
  std::unique_ptr<DistanceMatrixImpl> diagonal_impl_ptr_;    // PIMPL

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <class T> using NodeVector = std::vector<std::shared_ptr<PhyloNode<T>>>;

template<class T> class UPGMAMatrix {

public:

  explicit UPGMAMatrix(const NodeVector<T>& node_vector) : node_vector_(node_vector),
                                                           distance_matrix_(node_vector.size()) {
    initializeDistance();

  }
  ~UPGMAMatrix() = default;

  void calculateReduce();
  bool writeNewick(const std::string& file_name) const;

private:

  DistanceType_t distance(std::shared_ptr<PhyloNode<T>> row_node, std::shared_ptr<PhyloNode<T>> column_node) const;
  void initializeDistance();
  bool reduce(size_t row, size_t column, DistanceType_t minimum);
  void writeNode(std::shared_ptr<PhyloNode<T>> node, std::ofstream& newick_file) const;


  NodeVector<T> node_vector_;
  DistanceMatrix distance_matrix_;


};


template<class T>
DistanceType_t UPGMAMatrix<T>::distance(std::shared_ptr<PhyloNode<T>> row_node,
                                        std::shared_ptr<PhyloNode<T>> column_node) const {

  return row_node->leaf()->distance(column_node->leaf());

}


template<class T> void UPGMAMatrix<T>::initializeDistance() {

  for (size_t row = 0; row < node_vector_.size(); ++row) {

    for (size_t column = 0; column < row; column++) {

      distance_matrix_.setDistance(row, column, distance(node_vector_[row], node_vector_[column]));

    }

  }

}


template<class T> void UPGMAMatrix<T>::calculateReduce() {

  while (node_vector_.size() > 1) {

    size_t row;
    size_t column;
    DistanceType_t minimum = distance_matrix_.reduceMinimum(row, column);

    reduce(row, column, minimum);

  } // while reduce.

}


template<class T> bool UPGMAMatrix<T>::reduce(size_t row, size_t column, DistanceType_t minimum) {

  NodeVector<T> temp_node_vector;
  std::shared_ptr<PhyloNode<T>> column_node = nullptr;
  std::shared_ptr<PhyloNode<T>> row_node = nullptr;

  for (size_t idx = 0; idx < node_vector_.size(); idx++) {

    if (not (idx == column or idx == row)) {

      temp_node_vector.push_back(node_vector_[idx]);

    } else if (idx == column) {

      column_node = node_vector_[idx];

    } else if (idx == row) {

      row_node = node_vector_[idx];

    }

  }

  if (not column_node or not row_node) {

    ExecEnv::log().error("Null pointer found, error calculating UPGMA, row: {}, column: {}, nodes: {}",
                         row, column, node_vector_.size());
    return false;

  }

  node_vector_ = temp_node_vector;

  std::shared_ptr<const T> merged_leaf = row_node->leaf()->merge(column_node->leaf());
  DistanceType_t node_distance = minimum / 2;
  std::shared_ptr<PhyloNode<T>> merged_node(std::make_shared<PhyloNode<T>>(merged_leaf, node_distance));
  merged_node->addNode(row_node);
  merged_node->addNode(column_node);

  node_vector_.push_back(merged_node);

  return true;

}


template<class T> bool UPGMAMatrix<T>::writeNewick(const std::string& file_name) const {

  std::ofstream newick_file;

  // Open input file.

  newick_file.open(file_name);

  if (not newick_file.good()) {

    ExecEnv::log().error("I/O error; could not open Newick file: {}", file_name);
    return false;

  }

  for (auto node : node_vector_) {

    writeNode(node, newick_file);

  }

  newick_file << ";";

  newick_file.close();

  return true;

}


template<class T> void UPGMAMatrix<T>::writeNode(std::shared_ptr<PhyloNode<T>> node, std::ofstream& newick_file) const {


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

  }

  node->leaf()->write_node(newick_file);
  newick_file << ":";
  newick_file << node->distance();

}


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_STATISTICS_PHYLO_H
