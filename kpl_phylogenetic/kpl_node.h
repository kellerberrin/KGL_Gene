//
// Created by kellerberrin on 10/12/19.
//

#ifndef KPL_NODE_H
#define KPL_NODE_H

#include "kpl_splittree.h"

#include <string>
#include <vector>
#include  <iostream>


namespace kellerberrin::phylogenetic {   //  organization::project level namespace

class TreeIO;

class Node {


public:

  friend TreeIO;

  using PtrNode = Node*;
  using ConstPtrNode = Node const *;
  using PtrVector = std::vector<PtrNode>;
  using ConstPtrVector = std::vector<ConstPtrNode>;
  using Vector = std::vector<Node>;

  Node() { clear(); }
  Node(const Node& node) = default;
  ~Node() = default;

  Node& operator=(const Node& node) = default;

// Access
  [[nodiscard]] const PtrNode& getParent() const { return _parent; }
  [[nodiscard]] const PtrNode& getLeftChild() const { return _left_child; }
  [[nodiscard]] const PtrNode& getRightSib() const { return _right_sib; }
  [[nodiscard]] long getNumber() const { return _number; }
  [[nodiscard]] std::string getName() const { return _name; }
  [[nodiscard]] const Split& getConstSplit() { return _split; }
  [[nodiscard]] double getEdgeLength() const { return _edge_length; }
  [[nodiscard]] bool checkValidNumber() const { return getNumber() != _NO_NUMBER; }  // true if valid.

  [[nodiscard]] static double smallestEdgeLength() { return _SMALLEST_EDGE_LENGTH; }
  [[nodiscard]] static PtrNode nullNode() { return nullptr; }
  [[nodiscard]] static bool isNullNode(ConstPtrNode node) { return node == nullNode(); }
  [[nodiscard]] static std::string printNode(Node::ConstPtrNode Node);
  // Modify

  [[nodiscard]] Split& getSplit() { return _split; }
  void setParent(const PtrNode& parent) { _parent = parent; }
  void setLeftChild(const PtrNode& left_child) { _left_child = left_child; }
  void setRightSib(const PtrNode& right_sib) { _right_sib = right_sib; }
  void setEdgeLength(double v);
  void setName(const std::string& name) { _name = name; }
  void setNumber(long number) { _number = number; }
  void inValidNumber() { _number = _NO_NUMBER; }
// Flags
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


  void clearPointers();


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
  long _number;
  std::string _name;
  double _edge_length;
  Split _split;
  unsigned _flags;

  static constexpr const double _SMALLEST_EDGE_LENGTH = 1.0e-12; // This appears to be a "magic number".
  static constexpr const double _NO_NUMBER = -1; // Used as an initialization flag, todo::replace ASAP.

};


} // end namespace


#endif // KPL_NODE_H
