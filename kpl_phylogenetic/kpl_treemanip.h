//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_TREEMANIP_H
#define KPL_TREEMANIP_H


#include "kpl_xstrom.h"
#include "kpl_tree.h"
#include "kpl_splittree.h"


#include <boost/range/adaptor/reversed.hpp>
#include <boost/format.hpp>


#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#include <regex>


namespace kellerberrin::phylogenetic {   //  organization level namespace


class TreeManip {

public:

  TreeManip();
  explicit TreeManip(std::shared_ptr<Tree> t);
  ~TreeManip() = default;

///////////////////////////////////////////////////////////////////////////////////
// Static functions that perform the actual manipulations.

  static void sscaleAllEdgeLengths(std::shared_ptr<Tree> tree, double scaler);
  static void srefreshPreorder(std::shared_ptr<Tree> tree);
  static void srefreshLevelorder(std::shared_ptr<Tree> tree);
  static void srenumberInternals(std::shared_ptr<Tree> tree);
  static void srerootAtNodeNumber(std::shared_ptr<Tree> tree, int node_number);
  static void srerootAtNode(std::shared_ptr<Tree> tree, Node::PtrNode prospective_root);
  static void sstoreSplits(std::shared_ptr<Tree> tree, std::set<Split> &splitset);
  static void sselectAllPartials(std::shared_ptr<Tree> tree);
  static void sdeselectAllPartials(std::shared_ptr<Tree> tree);
  static void sselectAllTMatrices(std::shared_ptr<Tree> tree);
  static void sdeselectAllTMatrices(std::shared_ptr<Tree> tree);
  static void selectPartialsHereToRoot(Node::PtrNode  a);
  static void sflipPartialsAndTMatrices(std::shared_ptr<Tree> tree);
  static void sLargetSimonSwap(std::shared_ptr<Tree> tree, Node::PtrNode  a, Node::PtrNode  b);
  static void insertSubtreeOnLeft(Node::PtrNode  subtree, Node::PtrNode  parent);
  static void insertSubtreeOnRight(Node::PtrNode  subtree, Node::PtrNode  parent);
  static void detachSubtree(Node::PtrNode  subtree);
  static void srectifyNumInternals(std::shared_ptr<Tree> tree, int incr);
  static void srefreshNavigationPointers(std::shared_ptr<Tree> tree);
  static void sputUnusedNode(std::shared_ptr<Tree> tree, Node::PtrNode  node);

  [[nodiscard]] static double scalcTreeLength(std::shared_ptr<const Tree> tree);
  [[nodiscard]] static Node::PtrNode findNextPreorder(Node::ConstPtrNode  node);
  [[nodiscard]] static std::shared_ptr<Tree> screateTestTree();
  [[nodiscard]] static Node::PtrNode  srandomInternalEdge(std::shared_ptr<Tree> tree, double uniform_deviate);
  [[nodiscard]] static unsigned scountEdges(std::shared_ptr<const Tree> tree);
  [[nodiscard]] static unsigned scalcResolutionClass(std::shared_ptr<const Tree> tree);
  [[nodiscard]] static unsigned countChildren(Node::ConstPtrNode  nd);
  [[nodiscard]] static unsigned scountInternals(std::shared_ptr<const Tree> tree);
  [[nodiscard]] static Node::PtrNode findLeftSib(Node::ConstPtrNode  node);
  [[nodiscard]] static Node::PtrNode findRightmostChild(Node::ConstPtrNode  pNode);
  [[nodiscard]] static Node::PtrNode findLastPreorderInClade(Node::PtrNode  start);
  [[nodiscard]] static Node::PtrNode sgetUnusedNode(std::shared_ptr<Tree> tree);

////////////////////////////////////////////////////////////////////////////////////

  [[nodiscard]] Node::ConstPtrNode getConstNode(size_t node_index) const { return  getConstTree()->getConstNode(node_index); }
  [[nodiscard]] const std::shared_ptr<const Tree> getConstTree() const { return _tree; }

  [[nodiscard]] std::shared_ptr<Tree> getTree() { return _tree; }
  [[nodiscard]] Node::PtrNode getNode(size_t node_index) { return  getTree()->getNode(node_index); }

  [[nodiscard]] double calcTreeLength() const;
  [[nodiscard]] unsigned calcResolutionClass() const;
  [[nodiscard]] unsigned countEdges() const;
  [[nodiscard]] unsigned countInternals() const;
  [[nodiscard]] Node::PtrNode randomInternalEdge(double uniform01);
  [[nodiscard]] std::string makeNewick(unsigned precision) const;
  [[nodiscard]] Node::PtrNode getUnusedNode();

  void scaleAllEdgeLengths(double scaler);
  void createTestTree();
  void setTree(std::shared_ptr<Tree> tree) { assert(tree); _tree = tree; }
  void buildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies);
  void storeSplits(std::set<Split> & splitset);
  void rerootAtNodeNumber(int node_number);
  void LargetSimonSwap(Node::PtrNode  a, Node::PtrNode  b);
  void nniNodeSwap(Node::PtrNode  a, Node::PtrNode  b);
  void rectifyNumInternals(int incr);
  void refreshNavigationPointers();
  void putUnusedNode(Node::PtrNode  node);
  void selectAllPartials();
  void deselectAllPartials();
  void selectAllTMatrices();
  void deselectAllTMatrices();
  void flipPartialsAndTMatrices();
  void clear();

private:

  void refreshPreorder();
  void refreshLevelorder();
  void renumberInternals();
  void rerootAtNode(Node::PtrNode  prospective_root);

  std::shared_ptr<Tree> _tree;

public:

  typedef std::shared_ptr< TreeManip > SharedPtr;
};


} // end namespace


#endif // KPL_TREEMANIP_H
