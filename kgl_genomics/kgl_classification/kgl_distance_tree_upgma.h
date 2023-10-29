//
// Created by kellerberrin on 16/12/17.
//

#ifndef KGL_STATISTICS_UPGMA_H
#define KGL_STATISTICS_UPGMA_H

#include "kgl_distance_tree_base.h"
#include "kgl_runtime_resource.h"
#include "kel_utility.h"
#include "kgl_sequence_distance.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern to isolate Boost functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DistanceMatrix {

public:

  explicit DistanceMatrix();
  explicit DistanceMatrix(size_t matrix_size);
  explicit DistanceMatrix(const DistanceMatrix& copy);
  virtual ~DistanceMatrix();  // Do not use the default destructor, see PIMPL fwd decl below.


  [[nodiscard]] DistanceType_t minimum(size_t& i, size_t& j) const;
  [[nodiscard]] DistanceType_t maximum(size_t& i, size_t& j) const;
  void setDistance(size_t i, size_t j, DistanceType_t distance);
  [[nodiscard]] DistanceType_t getDistance(size_t i, size_t j) const;

  [[nodiscard]] size_t size() const;
  void resize(size_t new_size);

private:

  class DistanceMatrixImpl;       // Forward declaration of the boost strict diagonal implementation class
  std::unique_ptr<DistanceMatrixImpl> diagonal_impl_ptr_;    // PIMPL


};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DistanceTreeUPGMA : public DistanceTreeBase {

public:

  explicit DistanceTreeUPGMA() = default;
  ~DistanceTreeUPGMA() override = default;


  void calculateTree(std::shared_ptr<DistanceNodeVector> node_vector_ptr) override {

    node_vector_ptr_ = node_vector_ptr;
    distance_matrix_.resize(node_vector_ptr->size());
    initializeDistance();
    UPGMATree();

  }

  [[nodiscard]] bool writeNewick(const std::string& file_name) const override ;

private:

  [[nodiscard]] DistanceType_t distance(std::shared_ptr<TreeDistanceNode> row_node, std::shared_ptr<TreeDistanceNode> column_node) const;
  void initializeDistance();
  virtual void normalizeDistance();
  void rescaleDistance();
  void identityZeroDistance();
  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void reduceDistance(size_t i, size_t j);
  void writeNode(const std::shared_ptr<TreeDistanceNode>& node, std::ofstream& newick_file) const;
  [[nodiscard]] size_t getLeafCount(size_t leaf_idx) const;
  void UPGMATree();

  DistanceMatrix distance_matrix_;
  std::shared_ptr<DistanceNodeVector> node_vector_ptr_;

};


}   // end namespace


#endif //KGL_STATISTICS_UPGMA_H
