//
// Created by kellerberrin on 23/05/19.
//

#ifndef KGL_PHYLOGENETIC_TREE_H
#define KGL_PHYLOGENETIC_TREE_H


#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kel_exec_env.h"


namespace kellerberrin::genome {   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Virtual Distance class implemented elsewhere that actually calculates the Phylo tree calculateDistance.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using DistanceType_t = double;

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
  // Two identical distance node objects (same const void*) have zero calculateDistance.
  // This function is only re-defined and used if the calculateDistance metric needs to set a particular
  // condition for a zero distance. Most calculateDistance metrics will not need to re-define this function.
  [[nodiscard]] virtual bool zeroDistance(std::shared_ptr<const VirtualDistanceNode> node) const;

private:

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tree Classification Nodes.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
class PhyloNode;  // fwd.
using PhyloNodeVector = std::vector<std::shared_ptr<PhyloNode>>;
using OutNodes = std::multimap<DistanceType_t , std::shared_ptr<PhyloNode>>;

class PhyloNode {

public:

  explicit PhyloNode(std::shared_ptr<const VirtualDistanceNode> node) : node_(node), distance_(0.0) {}
  ~PhyloNode() = default;

  void addOutNode(std::shared_ptr<PhyloNode> node) {
    out_nodes_.insert(std::pair<DistanceType_t , std::shared_ptr<PhyloNode>>(node->distance(), node));
  }

  [[nodiscard]] DistanceType_t distance() const { return distance_; }
  void distance(DistanceType_t update) { distance_ = update; }

  [[nodiscard]] std::shared_ptr<const VirtualDistanceNode> node() const { return node_; }
  [[nodiscard]] const OutNodes& outNodes() const { return out_nodes_; }

  [[nodiscard]] bool isLeaf() const { return outNodes().empty(); }
  // Recursively counts the total number of leaf nodes, returns 1 if this is a node
  [[nodiscard]] size_t leafNodeCount() const;

private:

  std::shared_ptr<const VirtualDistanceNode> node_;
  DistanceType_t distance_;
  OutNodes out_nodes_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract Phylogenetic Tree
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class DistanceTree {

public:

  DistanceTree() {}
  virtual ~DistanceTree() = default;


  virtual void calculateTree(std::shared_ptr<PhyloNodeVector> node_vector_ptr) = 0;

  virtual bool writeNewick(const std::string &file_name) const = 0;

protected:


};




}   // end namespace


#endif //KGL_PHYLOGENETIC_TREE_H
