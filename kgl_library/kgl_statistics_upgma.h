//
// Created by kellerberrin on 16/12/17.
//

#ifndef KGL_STATISTICS_UPGMA_H
#define KGL_STATISTICS_UPGMA_H

#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kgl_exec_env.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db_population.h"
#include "kgl_utility.h"
#include "kgl_sequence_distance.h"
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual Distance class implemented elsewhere that actually calculates the UPGMA distance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DistanceType_t = double;
class PhyloNode;  // fwd.
using PhyloNodeVector = std::vector<std::shared_ptr<PhyloNode>>;

class UPGMADistanceNode {

public:

  UPGMADistanceNode() = default;
  UPGMADistanceNode(const UPGMADistanceNode&) = default;
  virtual ~UPGMADistanceNode() = default;

  // UPGMA Classification functions
  // Function to tag the nodes. Override as necessary.
  virtual void writeNode(std::ostream& outfile) const = 0;
  // Pure Virtual calculates the distance between nodes.
  virtual DistanceType_t distance(std::shared_ptr<const UPGMADistanceNode> distance_node) const = 0;
  // Pure Virtual calculates the zero distance between nodes.
  // This function is only re-defined and used if the distance metric needs to set a particular
  // condition for a zero distance. Most distance metrics will not need to re-define this function.
  virtual DistanceType_t zeroDistance(std::shared_ptr<const UPGMADistanceNode>) const { return 1.0; }

private:

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<PhyloNode>>;
class PhyloNode {

public:

  explicit PhyloNode(std::shared_ptr<const UPGMADistanceNode> leaf) : leaf_(leaf), distance_(0) {}
  ~PhyloNode() = default;

  void addOutNode(std::shared_ptr<PhyloNode> node) {
    out_nodes_.insert(std::pair<DistanceType_t , std::shared_ptr<PhyloNode>>(node->distance(), node));
  }

  DistanceType_t distance() const { return distance_; }
  void distance(DistanceType_t update) { distance_ = update; }

  std::shared_ptr<const UPGMADistanceNode> leaf() const { return leaf_; }
  const OutNodes& getMap() const { return out_nodes_; }

  // Recursively counts the total number of leaf nodes.
  bool isLeaf() const { return getMap().empty(); }
  size_t leafNodeCount() const;

private:

  std::shared_ptr<const UPGMADistanceNode> leaf_;
  DistanceType_t distance_;
  OutNodes out_nodes_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Distance matrix. Implements the PIMPL pattern to isolate Boost functionality.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DistanceMatrix {

public:

  explicit DistanceMatrix(size_t matrix_size);
  explicit DistanceMatrix(const DistanceMatrix& copy);
  virtual ~DistanceMatrix();  // Do not use the default destructor, see PIMPL fwd decl below.


  DistanceType_t minimum(size_t& i, size_t& j) const;
  DistanceType_t maximum(size_t& i, size_t& j) const;
  void reduce(size_t i, size_t j);
  void setDistance(size_t i, size_t j, DistanceType_t distance);
  DistanceType_t getDistance(size_t i, size_t j) const;

  virtual size_t getLeafCount(size_t) const { return 1; }

private:

  class BoostDistanceMatrix;       // Forward declaration of the boost strict diagonal implementation class
  std::unique_ptr<BoostDistanceMatrix> diagonal_impl_ptr_;    // PIMPL

  size_t size() const;
  void resize(size_t new_size);

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UPGMA Distance matrix
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//using MatrixNode = PhyloNode;

class UPGMAMatrix : DistanceMatrix {

public:

  explicit UPGMAMatrix(std::shared_ptr<PhyloNodeVector> node_vector_ptr) : DistanceMatrix(node_vector_ptr->size()),
                                                                           node_vector_ptr_(node_vector_ptr) {
    initializeDistance();

  }
  ~UPGMAMatrix() override = default;


  void calculateReduce();

  bool writeNewick(const std::string& file_name) const;

private:

  DistanceType_t distance(std::shared_ptr<PhyloNode> row_node, std::shared_ptr<PhyloNode> column_node) const;
  DistanceType_t zeroDistance(std::shared_ptr<PhyloNode> row_node, std::shared_ptr<PhyloNode> column_node) const;
  void initializeDistance();
  virtual void normalizeDistance();
  void rescaleDistance();
  void identityZeroDistance();
  bool reduceNode(size_t row, size_t column, DistanceType_t minimum);
  void writeNode(std::shared_ptr<PhyloNode> node, std::ofstream& newick_file) const;
  size_t getLeafCount(size_t leaf_idx) const override;

  std::shared_ptr<PhyloNodeVector> node_vector_ptr_;


};


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_STATISTICS_UPGMA_H
