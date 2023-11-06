//
// Created by kellerberrin on 23/05/19.
//

#ifndef KGL_DISTANCE_TREE_H
#define KGL_DISTANCE_TREE_H


#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kel_exec_env.h"
#include "kgl_distance_matrix.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual Distance class implemented elsewhere that actually calculates the Phylo tree calculateDistance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class VirtualDistanceNode {

public:

  VirtualDistanceNode() = default;
  VirtualDistanceNode(const VirtualDistanceNode&) = default;
  virtual ~VirtualDistanceNode() = default;

  // Classification functions
  // Function to tag the nodes. Override as necessary.
  virtual void writeNode(std::ostream& outfile) const = 0;
  // Pure Virtual calculates the calculateDistance between nodes.
  [[nodiscard]] virtual DistanceType_t distance(std::shared_ptr<const VirtualDistanceNode> distance_node) const = 0;

private:

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
class TreeDistanceNode;  // fwd.
using DistanceNodeVector = std::vector<std::shared_ptr<TreeDistanceNode>>;
using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<TreeDistanceNode>>;

class TreeDistanceNode {

public:

  explicit TreeDistanceNode(std::shared_ptr<const VirtualDistanceNode> node) : node_(std::move(node)), distance_(0.0) {}
  ~TreeDistanceNode() = default;

  void addOutNode(const std::shared_ptr<TreeDistanceNode>& node) { out_nodes_.insert({node->distance(), node}); }

  [[nodiscard]] DistanceType_t distance() const { return distance_; }
  void distance(DistanceType_t update) { distance_ = update; }

  [[nodiscard]] std::shared_ptr<const VirtualDistanceNode> node() const { return node_; }
  [[nodiscard]] const OutNodes& outNodes() const { return out_nodes_; }

  [[nodiscard]] bool isLeaf() const { return outNodes().empty(); }
  // Recursively counts the total number of leaf nodes, returns 1 if this is a leaf node
  [[nodiscard]] size_t leafNodeCount() const;

private:

  std::shared_ptr<const VirtualDistanceNode> node_;
  DistanceType_t distance_;
  OutNodes out_nodes_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract Phylogenetic Tree
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceTreeBase {

public:

  DistanceTreeBase() = default;
  virtual ~DistanceTreeBase() = default;


  virtual void calculateTree(std::shared_ptr<DistanceNodeVector> node_vector_ptr) = 0;

  [[nodiscard]] virtual bool writeNewick(const std::string &file_name) const = 0;

private:


};




}   // end namespace


#endif //KGL_DISTANCE_TREE_H
