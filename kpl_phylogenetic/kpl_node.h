//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_NODE_H
#define KPL_NODE_H

#include "kpl_splittree.h"

#include <string>
#include <vector>
#include  <iostream>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


//class Tree;
//class TreeManip;
class Likelihood;
//class Updater;

class Node {

  friend class Tree;

  friend class TreeManip;
  friend class Likelihood;
//  friend class Updater;

public:

  using PtrNode = Node*;
  using PtrVector = std::vector<PtrNode>;
  using Vector = std::vector<Node>;

  Node();
  ~Node();

  PtrNode getParent() { return _parent; }
  PtrNode getLeftChild() { return _left_child; }
  PtrNode getRightSib() { return _right_sib; }

  int getNumber() { return _number; }

  std::string getName() { return _name; }

  Split getSplit() { return _split; }

  bool                isSelected() { return _flags & Flag::Selected; }
  void                select() { _flags |= Flag::Selected; }
  void                deselect() { _flags &= ~Flag::Selected; }


  bool                isSelPartial() { return _flags & Flag::SelPartial; }
  void                selectPartial() { _flags |= Flag::SelPartial; }
  void                deselectPartial() { _flags &= ~Flag::SelPartial; }


  bool                isSelTMatrix() { return _flags & Flag::SelTMatrix; }
  void                selectTMatrix() { _flags |= Flag::SelTMatrix; }
  void                deselectTMatrix() { _flags &= ~Flag::SelTMatrix; }


  bool                isAltPartial() { return _flags & Flag::AltPartial; }
  void                setAltPartial() { _flags |= Flag::AltPartial; }
  void                clearAltPartial() { _flags &= ~Flag::AltPartial; }


  bool                isAltTMatrix() { return _flags & Flag::AltTMatrix; }
  void                setAltTMatrix() { _flags |= Flag::AltTMatrix; }
  void                clearAltTMatrix() { _flags &= ~Flag::AltTMatrix; }

  [[nodiscard]] double getEdgeLength() const { return _edge_length; }

  void setEdgeLength(double v);

  void clearPointers() { _left_child = _right_sib = _parent = nullptr; }

  static const double _smallest_edge_length;

private:

  enum Flag {
    Selected   = (1 << 0),
    SelPartial = (1 << 1),
    SelTMatrix = (1 << 2),
    AltPartial = (1 << 3),
    AltTMatrix = (1 << 4)
  };

  void clear();

  PtrNode _left_child;
  PtrNode _right_sib;
  PtrNode _parent;
  int _number;
  std::string _name;
  double _edge_length;
  Split _split;
  unsigned _flags;

};


} // phylogenetic
} // kellerberrin


#endif // KPL_NODE_H
