//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_TREE_H
#define KPL_TREE_H

#include "kpl_node.h"

#include <memory>
#include <iostream>
#include <stack>

namespace kellerberrin::phylogenetic {   //  organization level namespace

class TreeManip;

// Only TreeManip can modify this object.
class Tree {

public:

  friend TreeManip;

  explicit Tree(unsigned num_nodes = 0);
  ~Tree() = default;


  // Access
  [[nodiscard]] bool isRooted() const { return _rooted; } // Flag to indicate if the tree is rooted.
  [[nodiscard]] unsigned numLeaves() const { return _nleaves; }
  [[nodiscard]] unsigned numInternals() const { return _ninternals; }
  [[nodiscard]] unsigned numNodes() const { return static_cast<unsigned>(_nodes.size());}
  [[nodiscard]] const Node::PtrVector& getConstPreOrder() const { return _preorder; }
  [[nodiscard]] const Node::PtrVector& getConstLevelOrder() const { return _levelorder; }
  [[nodiscard]] const Node::Vector& getConstNodes() const { return _nodes; }
  [[nodiscard]] const std::stack<Node::PtrNode >& getUnUsed() const { return _unused_nodes; }
  [[nodiscard]] Node::ConstPtrNode getConstNode(size_t node_index) const;
  // The top (rooted) node of the tree. Does NOT indicate the tree is rooted.
  [[nodiscard]] Node::ConstPtrNode getConstRoot() const { return _root; }
  [[nodiscard]] std::string treeDescription() const;

  // Modify
  [[nodiscard]] Node::PtrNode getNode(size_t node_index);    // Only used in TreeManip.
  [[nodiscard]] Node::PtrNode getRootNode() const { return _root; }
  [[nodiscard]] Node::PtrVector getNodes();
  [[nodiscard]] const Node::PtrVector& getPreOrder() const { return _preorder; }
  void clearPreOrder() { _preorder.clear(); }
  void pushPreOrder(Node::PtrNode node) { _preorder.push_back(node); }
  void clearLevelOrder() { _levelorder.clear(); }
  void pushLevelOrder(Node::PtrNode node) { _levelorder.push_back(node); }
  void pushUnused(Node::PtrNode node) { _unused_nodes.push(node); }
  void popUnused() { _unused_nodes.pop(); }
  void setRooted(bool rooted) { _rooted = rooted; }
  // The top (rooted) node of the tree. Does NOT indicate the tree is rooted.
  void setRootNode(Node::PtrNode root) { _root = root; }
  void setLeaves(unsigned leaves) { _nleaves = leaves; }
  void setInternals(unsigned internals) { _ninternals = internals; }


private:

  void clear();

  bool _rooted;
  Node::PtrNode _root;
  unsigned _nleaves;
  unsigned _ninternals;
  Node::PtrVector _preorder;
  Node::PtrVector _levelorder;
  Node::Vector _nodes;
  std::stack<Node::PtrNode >  _unused_nodes;


};



}   // end namespace


#endif // KPL_TREE_H
