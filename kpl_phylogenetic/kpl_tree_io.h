//
// Created by kellerberrin on 12/1/20.
//

#ifndef KPL_TREE_IO_H
#define KPL_TREE_IO_H

#include "kpl_tree.h"

#include <memory>
#include <iostream>
#include <stack>

namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


// Performs TreeIO to and from the Newick tree format.

class TreeIO {

public:

  TreeIO() = delete;
  ~TreeIO() = delete;

  [[nodiscard]] static Tree::SharedPtr sbuildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies);
  [[nodiscard]] static std::string smakeNewick(Tree::ConstSharedPtr tree, unsigned precision);
  [[nodiscard]] static Tree::SharedPtr buildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies);

private:

  static void extractNodeNumberFromName(Node::PtrNode  node, std::set<unsigned> & used);
  static void extractEdgeLen(Node::PtrNode nd, std::string edge_length_string);
  static unsigned countNewickLeaves(const std::string& newick);
  static void stripOutNexusComments(std::string& newick);
  static bool canHaveSibling(Node::PtrNode node, bool rooted, bool allow_polytomies);


};


}
}

#endif //KGL_KPL_TREE_IO_H
