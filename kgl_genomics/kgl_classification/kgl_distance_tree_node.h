//
// Created by kellerberrin on 28/12/23.
//

#ifndef KGL_DISTANCE_TREE_NODE_H
#define KGL_DISTANCE_TREE_NODE_H



#include <memory>
#include <map>
#include <vector>
#include <fstream>

#include "kel_exec_env.h"
#include "kgl_distance_matrix.h"


namespace kellerberrin::genome {   //  organization level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Tree Classification Nodes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Used by the classification functions.
class TreeNodeDistance;  // fwd.
using BaseOutNodes = std::multimap<DistanceType_t , std::shared_ptr<TreeNodeDistance>>;
using BaseParentNode = std::shared_ptr<TreeNodeDistance>;

class TreeNodeDistance {

public:

  TreeNodeDistance() = default;
  virtual ~TreeNodeDistance() = default;

  void addOutNode(const std::shared_ptr<TreeNodeDistance>& node) { out_nodes_.insert({node->parentDistance(), node}); }

  [[nodiscard]] DistanceType_t parentDistance() const { return parent_distance_; }
  void parentDistance(DistanceType_t update) { parent_distance_ = update; }
  [[nodiscard]] const std::optional<BaseParentNode>& parentNode() const { return parent_node_; }
  void parentNode(const std::shared_ptr<TreeNodeDistance>& parent_ptr) { parent_node_ = parent_ptr; }

  [[nodiscard]] const BaseOutNodes& outNodes() const { return out_nodes_; }

  // Node type.
  [[nodiscard]] bool isRoot() const { return not parent_node_; }
  [[nodiscard]] bool isLeaf() const { return outNodes().empty(); }
  [[nodiscard]] bool isBranch() const { return not isRoot() and not isLeaf(); }
  // By default, recursively counts the total number of leaf nodes, returns 1 if this is a leaf node
  [[nodiscard]] virtual size_t leafNodeCount() const;
  // Virtual text and distance.
  [[nodiscard]] virtual std::string nodeText() const = 0;
  [[nodiscard]] virtual DistanceType_t distance(const std::shared_ptr<const TreeNodeDistance>& distance_node) const;

private:

  BaseOutNodes out_nodes_;
  DistanceType_t parent_distance_{0.0};
  std::optional<BaseParentNode> parent_node_{std::nullopt};

};

using TreeNodeVector = std::vector<std::shared_ptr<TreeNodeDistance>>;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Clade (Branch) Nodes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CladeNode : public TreeNodeDistance {

public:

  CladeNode(std::string clade_description, size_t leaf_nodes)
      : clade_description_(std::move(clade_description)), leaf_nodes_(leaf_nodes) {}
  ~CladeNode() override = default;

  // Classification functions
  [[nodiscard]] std::string nodeText() const override { return clade_description_; }
  [[nodiscard]] size_t leafNodeCount() const override { return leaf_nodes_; }

private:

  std::string clade_description_;
  size_t leaf_nodes_;

};



}   // end namespace


#endif //KGL_DISTANCE_TREE_NODE_H
