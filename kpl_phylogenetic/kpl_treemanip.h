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


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class TreeManip {

public:

  TreeManip();
  TreeManip(Tree::SharedPtr t);
  ~TreeManip();

///////////////////////////////////////////////////////////////////////////////////
// Static functions that perform the actual manipulations.
  static double scalcTreeLength(Tree::ConstSharedPtr tree);
  static void sscaleAllEdgeLengths(Tree::SharedPtr tree, double scaler);
  [[nodiscard]] static Tree::SharedPtr screateTestTree();
  static void srefreshPreorder(Tree::SharedPtr tree);
  [[nodiscard]] static Node::PtrNode findNextPreorder(Node::ConstPtrNode  node);
  static void srefreshLevelorder(Tree::SharedPtr tree);
  static void srenumberInternals(Tree::SharedPtr tree);
  static void srerootAtNodeNumber(Tree::SharedPtr tree, int node_number);
  static void srerootAtNode(Tree::SharedPtr tree, Node::PtrNode prospective_root);
  static void sstoreSplits(Tree::SharedPtr tree, std::set<Split> &splitset);
  static void sselectAllPartials(Tree::SharedPtr tree);
  static void sdeselectAllPartials(Tree::SharedPtr tree);
  static void sselectAllTMatrices(Tree::SharedPtr tree);
  static void sdeselectAllTMatrices(Tree::SharedPtr tree);
  static void selectPartialsHereToRoot(Node::PtrNode  a);
  static void sflipPartialsAndTMatrices(Tree::SharedPtr tree);
  static void sLargetSimonSwap(Tree::SharedPtr tree, Node::PtrNode  a, Node::PtrNode  b);
  [[nodiscard]] static Node::PtrNode  srandomInternalEdge(Tree::SharedPtr tree, double uniform_deviate);
  [[nodiscard]] static unsigned scountEdges(Tree::ConstSharedPtr tree);
  [[nodiscard]] static unsigned scalcResolutionClass(Tree::ConstSharedPtr tree);
  [[nodiscard]] static unsigned countChildren(Node::ConstPtrNode  nd);
  [[nodiscard]] static unsigned scountInternals(Tree::ConstSharedPtr tree);
  [[nodiscard]] static Node::PtrNode findLeftSib(Node::ConstPtrNode  node);
  [[nodiscard]] static Node::PtrNode findRightmostChild(Node::ConstPtrNode  node);
  [[nodiscard]] static Node::PtrNode findLastPreorderInClade(Node::PtrNode  start);
  static void insertSubtreeOnLeft(Node::PtrNode  subtree, Node::PtrNode  parent);
  static void insertSubtreeOnRight(Node::PtrNode  subtree, Node::PtrNode  parent);
  static void detachSubtree(Node::PtrNode  subtree);
  static void srectifyNumInternals(Tree::SharedPtr tree, int incr);
  static void srefreshNavigationPointers(Tree::SharedPtr tree);
  [[nodiscard]] static Node::PtrNode sgetUnusedNode(Tree::SharedPtr tree);
  static void sputUnusedNode(Tree::SharedPtr tree, Node::PtrNode  node);


////////////////////////////////////////////////////////////////////////////////////

  [[nodiscard]] Node::ConstPtrNode getConstNode(size_t node_index) const { return  getConstTree()->getConstNode(node_index); }
  [[nodiscard]] Tree::ConstSharedPtr getConstTree() const { return _tree; }

  void setTree(Tree::SharedPtr tree) { assert(tree); _tree = tree; }
  [[nodiscard]] Tree::SharedPtr getTree() { return _tree; }
  [[nodiscard]] Node::PtrNode getNode(size_t node_index) { return  getTree()->getNode(node_index); }

  double calcTreeLength() const;
  unsigned calcResolutionClass() const;
  unsigned countEdges() const;
  unsigned countInternals() const;
  void scaleAllEdgeLengths(double scaler);

  void createTestTree();
  std::string makeNewick(unsigned precision) const;

  void                        buildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies);
  void                        storeSplits(std::set<Split> & splitset);
  void                        rerootAtNodeNumber(int node_number);

  void                        LargetSimonSwap(Node::PtrNode  a, Node::PtrNode  b);
  Node::PtrNode                       randomInternalEdge(double uniform01);

  void                        nniNodeSwap(Node::PtrNode  a, Node::PtrNode  b);
  void                        rectifyNumInternals(int incr);
  void                        refreshNavigationPointers();
  Node::PtrNode                       getUnusedNode();
  void                        putUnusedNode(Node::PtrNode  node);

  void                        selectAllPartials();
  void                        deselectAllPartials();
  void                        selectAllTMatrices();
  void                        deselectAllTMatrices();

  void                        flipPartialsAndTMatrices();

  void                        clear();

private:

  void                        refreshPreorder();
  void                        refreshLevelorder();
  void                        renumberInternals();
  void                        rerootAtNode(Node::PtrNode  prospective_root);

  Tree::SharedPtr             _tree;

public:

  typedef std::shared_ptr< TreeManip > SharedPtr;
};




} // phylogenetic
} // kellerberrin


#endif // KPL_TREEMANIP_H
