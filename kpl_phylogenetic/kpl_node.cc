//
// Created by kellerberrin on 12/12/19.
//


#include "kpl_node.h"


namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here

kpl::Node::Node() {

//  std::cout << "Creating Node object" << std::endl;
  clear();

}


kpl::Node::~Node() {

//  std::cout << "Destroying Node object" << std::endl;

}


void kpl::Node::clear() {

  _flags = 0;
  clearPointers();
  _number = -1;
  _name = "";
  _edge_length = _smallest_edge_length;

}


void kpl::Node::setEdgeLength(double v) {

  _edge_length = (v < _smallest_edge_length ? _smallest_edge_length : v);

}

