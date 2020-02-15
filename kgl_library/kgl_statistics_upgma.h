//
// Created by kellerberrin on 16/12/17.
//

#ifndef KGL_STATISTICS_UPGMA_H
#define KGL_STATISTICS_UPGMA_H

#include "kgl_phylogenetic_tree.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_population.h"
#include "kgl_utility.h"
#include "kgl_sequence_distance.h"
#include "kgl_variant_db.h"


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


  DistanceType_t minimum(size_t& i, size_t& j) const;
  DistanceType_t maximum(size_t& i, size_t& j) const;
  void setDistance(size_t i, size_t j, DistanceType_t distance);
  DistanceType_t getDistance(size_t i, size_t j) const;

  size_t size() const;
  void resize(size_t new_size);

private:

  class BoostDistanceMatrix;       // Forward declaration of the boost strict diagonal implementation class
  std::unique_ptr<BoostDistanceMatrix> diagonal_impl_ptr_;    // PIMPL


};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class UPGMAMatrix : public DistanceTree {

public:

  explicit UPGMAMatrix() = default;
  ~UPGMAMatrix() override = default;


  void calculateTree(std::shared_ptr<PhyloNodeVector> node_vector_ptr) override {

    node_vector_ptr_ = node_vector_ptr;
    distance_matrix_.resize(node_vector_ptr->size());
    initializeDistance();
    UPGMATree();

  }

  bool writeNewick(const std::string& file_name) const override ;

private:

  DistanceType_t distance(std::shared_ptr<PhyloNode> row_node, std::shared_ptr<PhyloNode> column_node) const;
  void initializeDistance();
  virtual void normalizeDistance();
  void rescaleDistance();
  void identityZeroDistance();
  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void reduceDistance(size_t i, size_t j);
  void writeNode(const std::shared_ptr<PhyloNode>& node, std::ofstream& newick_file) const;
  size_t getLeafCount(size_t leaf_idx) const;
  void UPGMATree();

  DistanceMatrix distance_matrix_;
  std::shared_ptr<PhyloNodeVector> node_vector_ptr_;

};


}   // end namespace


#endif //KGL_STATISTICS_UPGMA_H
