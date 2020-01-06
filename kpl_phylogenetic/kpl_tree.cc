//
// Created by kellerberrin on 12/12/19.
//

#include "kpl_tree.h"

namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here

kpl::Tree::Tree() {

//  std::cout << "Constructing a Tree" << std::endl;
  clear();

}


kpl::Tree::~Tree() {
//  std::cout << "Destroying a Tree" << std::endl;
}


void kpl::Tree::clear() {

  _is_rooted = false;
  _root = 0;
  _nodes.clear();
  _preorder.clear();
  _levelorder.clear();

}


bool kpl::Tree::isRooted() const {

  return _is_rooted;

}


unsigned kpl::Tree::numLeaves() const {

  return _nleaves;

}


unsigned kpl::Tree::numInternals() const {

  return _ninternals;

}


unsigned kpl::Tree::numNodes() const {

  return (unsigned) _nodes.size();

}
