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

  void                        setTree(Tree::SharedPtr t);
  Tree::SharedPtr             getTree();

  double                      calcTreeLength() const;
  unsigned                    calcResolutionClass() const;
  unsigned                    countEdges() const;
  unsigned                    countInternals() const;
  void                        scaleAllEdgeLengths(double scaler);

  void                        createTestTree();
  std::string                 makeNewick(unsigned precision) const;

  void                        buildFromNewick(const std::string newick, bool rooted, bool allow_polytomies);
  void                        storeSplits(std::set<Split> & splitset);
  void                        rerootAtNodeNumber(int node_number);

  void                        LargetSimonSwap(Node::PtrNode  a, Node::PtrNode  b);
  Node::PtrNode                       randomInternalEdge(double uniform01);

  void                        nniNodeSwap(Node::PtrNode  a, Node::PtrNode  b);
  unsigned                    countChildren(Node::PtrNode  nd) const;
  Node::PtrNode                       findLeftSib(Node::PtrNode  nd);
  Node::PtrNode                       findNextPreorder(Node::PtrNode  nd);
  Node::PtrNode                       findRightmostChild(Node::PtrNode  nd);
  Node::PtrNode                       findLastPreorderInClade(Node::PtrNode  start);
  void                        insertSubtreeOnLeft(Node::PtrNode  s, Node::PtrNode  u);
  void                        insertSubtreeOnRight(Node::PtrNode  s, Node::PtrNode  u);
  void                        detachSubtree(Node::PtrNode  s);
  void                        rectifyNumInternals(int incr);
  void                        refreshNavigationPointers();
  Node::PtrNode                       getUnusedNode();
  void                        putUnusedNode(Node::PtrNode  nd);

  void                        selectAll();
  void                        deselectAll();
  void                        selectAllPartials();
  void                        deselectAllPartials();
  void                        selectAllTMatrices();
  void                        deselectAllTMatrices();

  void                        selectPartialsHereToRoot(Node::PtrNode  a);
  void                        flipPartialsAndTMatrices();

  void                        clear();

private:

  void                        refreshPreorder();
  void                        refreshLevelorder();
  void                        renumberInternals();
  void                        rerootAtNode(Node::PtrNode  prospective_root);
  void                        extractNodeNumberFromName(Node::PtrNode  nd, std::set<unsigned> & used);
  void                        extractEdgeLen(Node::PtrNode  nd, std::string edge_length_string);
  unsigned                    countNewickLeaves(const std::string newick);
  void                        stripOutNexusComments(std::string & newick);
  bool                        canHaveSibling(Node::PtrNode  nd, bool rooted, bool allow_polytomies);

  Tree::SharedPtr             _tree;

public:

  typedef std::shared_ptr< TreeManip > SharedPtr;
};




} // phylogenetic
} // kellerberrin


#endif // KPL_TREEMANIP_H
