//
// Created by kellerberrin on 12/12/19.
//


#include "kpl_node.h"

#include <sstream>

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


std::string kpl::Node::printNode(ConstPtrNode Node) {

  std::stringstream ss;

  if (not Node::isNullNode(Node)) {

    ss << "| Address:" << static_cast<const void*>(Node) << " Num:" << Node->getNumber() << ", Name:" << Node->getName() << ", Length:" << Node ->getEdgeLength();

    ss << ", Left Child: ";

    if (not Node::isNullNode(Node->getLeftChild())) {

      ss << static_cast<const void*>(Node->getLeftChild());

    } else {

      ss << "NULL";

    }

    ss << ", Right Sib: ";

    if (not Node::isNullNode(Node->getRightSib())) {

      ss << static_cast<const void*>(Node->getRightSib());

    } else {

      ss << "NULL";

    }

    ss << ", Parent: ";

    if (not Node::isNullNode(Node->getParent())) {

      ss  << static_cast<const void*>(Node->getParent());

    } else {

      ss << "NULL";

    }

    ss << " |";

  } else {

    ss << "| NULL |";

  }

  return ss.str();

}

