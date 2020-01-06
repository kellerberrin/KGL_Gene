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


//class TreeManip;
class Likelihood;
//class Updater;
class PolytomyUpdater;



class Tree {

  friend class TreeManip;
  friend class Likelihood;
//  friend class Updater;
  friend class PolytomyUpdater;

public:

  Tree();

  ~Tree();

  bool isRooted() const;

  unsigned numLeaves() const;

  unsigned numInternals() const;

  unsigned numNodes() const;

  [[nodiscard]] const Node::PtrVector& preOrder() const { return _preorder; }

private:

  void clear();

  bool _is_rooted;
  Node::PtrNode _root;
  unsigned _nleaves;
  unsigned _ninternals;
  Node::PtrVector _preorder;
  Node::PtrVector _levelorder;
  Node::Vector _nodes;
  std::stack<Node::PtrNode >  _unused_nodes;

public:

  typedef std::shared_ptr<Tree> SharedPtr;

};



}   // phylogenetic
}  // kellerberrin


#endif // KPL_TREE_H
