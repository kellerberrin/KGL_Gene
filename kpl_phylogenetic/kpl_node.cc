//
// Created by kellerberrin on 12/12/19.
//


#include "kpl_node.h"


namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here


void kpl::Node::clear() {

  _flags = 0;
  inValidNumber();
  _name = "";
  _edge_length = _SMALLEST_EDGE_LENGTH;
  _left_child = nullNode();
  _right_sib = nullNode();
  _parent = nullNode();

}


void kpl::Node::clearPointers() {

  _left_child = nullNode();
  _right_sib = nullNode();
  _parent = nullNode();

}



void kpl::Node::setEdgeLength(double v) {

  _edge_length = (v < _SMALLEST_EDGE_LENGTH ? _SMALLEST_EDGE_LENGTH : v);

}


