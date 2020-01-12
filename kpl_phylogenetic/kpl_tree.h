//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_TREE_H
#define KPL_TREE_H

#include "kpl_node.h"

#include <memory>
#include <iostream>
#include <stack>

namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace




class Tree {

public:

  using SharedPtr = std::shared_ptr<Tree> ;
  using ConstSharedPtr = std::shared_ptr<const Tree> ;
  // Only TreeManip can modify this object.


  explicit Tree(unsigned num_nodes = 0);
  ~Tree() = default;

// Access
  [[nodiscard]] bool isRooted() const;
  [[nodiscard]] unsigned numLeaves() const;
  [[nodiscard]] unsigned numInternals() const;
  [[nodiscard]] unsigned numNodes() const;
  [[nodiscard]] const Node::PtrVector& getConstPreOrder() const { return _preorder; }
  [[nodiscard]] const Node::PtrVector& getConstLevelOrder() const { return _levelorder; }
  [[nodiscard]] Node::ConstPtrVector getConstNodes() const;
  [[nodiscard]] const std::stack<Node::PtrNode >& getUnUsed() const { return _unused_nodes; }
  [[nodiscard]] Node::ConstPtrNode getConstNode(size_t node_index) const;
  [[nodiscard]] Node::ConstPtrNode getConstRoot() const { return _root; }
// Modify
  [[nodiscard]] Node::PtrNode getNode(size_t node_index);    // Only used in TreeManip.
  [[nodiscard]] Node::PtrNode getRoot() { return _root; }
  [[nodiscard]] Node::PtrVector getNodes();
  [[nodiscard]] Node::PtrVector& getPreOrder() { return _preorder; }
  void clearPreOrder() { _preorder.clear(); }    // Only used in TreeManip.
  void pushPreOrder(Node::PtrNode node) { _preorder.push_back(node); }    // Only used in TreeManip.
  void clearLevelOrder() { _levelorder.clear(); }    // Only used in TreeManip.
  void pushLevelOrder(Node::PtrNode node) { _levelorder.push_back(node); }    // Only used in TreeManip.
  void pushUnused(Node::PtrNode node) { _unused_nodes.push(node); }    // Only used in TreeManip.
  void popUnused() { _unused_nodes.pop(); }    // Only used in TreeManip.
  void setRoot(Node::PtrNode root) { _root = root; } // Only used in TreeManip.
  void setLeaves(unsigned leaves) { _nleaves = leaves; } // Only used in TreeManip.
  void setInternals(unsigned internals) { _ninternals = internals; } // Only used in TreeManip.

private:

  void clear();

  Node::PtrNode _root;
  unsigned _nleaves;
  unsigned _ninternals;
  Node::PtrVector _preorder;
  Node::PtrVector _levelorder;
  Node::Vector _nodes;
  std::stack<Node::PtrNode >  _unused_nodes;


};



}   // phylogenetic
}  // kellerberrin


#endif // KPL_TREE_H
