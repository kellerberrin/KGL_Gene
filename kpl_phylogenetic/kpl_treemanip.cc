//
// Created by kellerberrin on 12/12/19.
//


#include "kpl_treemanip.h"
#include "kpl_tree_io.h"

#include  <cmath>

namespace kpl = kellerberrin::phylogenetic;

// This is where function bodies go

kpl::TreeManip::TreeManip() {
//  std::cerr << "Constructing a TreeManip" << std::endl;
  clear();
}


kpl::TreeManip::TreeManip(Tree::SharedPtr t) {
//  std::cerr << "Constructing a TreeManip with a supplied tree" << std::endl;
  clear();
  setTree(t);
}


kpl::TreeManip::~TreeManip() {
//  std::cerr << "Destroying a TreeManip" << std::endl;
}


void kpl::TreeManip::clear() {

  _tree.reset();

}


double kpl::TreeManip::calcTreeLength() const {

  return scalcTreeLength(getConstTree());

}

double kpl::TreeManip::scalcTreeLength(Tree::ConstSharedPtr tree) {

  double TL = 0.0;

  for (auto nd : tree->getConstPreOrder()) {

    TL += nd->getEdgeLength();

  }

  return TL;

}


void kpl::TreeManip::scaleAllEdgeLengths(double scaler) {

  sscaleAllEdgeLengths(getTree(), scaler);

}

void kpl::TreeManip::sscaleAllEdgeLengths(Tree::SharedPtr tree, double scaler) {

  for (auto nd : tree->getConstPreOrder()) {

    nd->setEdgeLength(scaler * nd->getEdgeLength());

  }

}


void kpl::TreeManip::createTestTree() {

  clear();
  setTree(screateTestTree());

}

kpl::Tree::SharedPtr kpl::TreeManip::screateTestTree() {

  Tree::SharedPtr tree(std::make_shared<Tree>(6));

  Node::PtrNode root_node = tree->getNode(0);
  Node::PtrNode first_internal = tree->getNode(1);
  Node::PtrNode second_internal = tree->getNode(2);
  Node::PtrNode first_leaf = tree->getNode(3);
  Node::PtrNode second_leaf = tree->getNode(4);
  Node::PtrNode third_leaf = tree->getNode(5);

  // Here is the structure of the tree (numbers in
  // parentheses are node numbers, other numbers
  // are edge lengths):
  //
  // first_leaf (0)   second_leaf (1)   third_leaf (2)
  //      \              /                  /
  //       \ 0.1        / 0.1              /
  //        \          /                  /
  //     second_internal (3)             / 0.2
  //             \                      /
  //              \ 0.1                /
  //               \                  /
  //                first_internal (4)
  //                        |
  //                        | 0.1
  //                        |
  //                    root_node (5)
  //
  root_node->setParent(Node::nullNode());
  root_node->setLeftChild(first_internal);
  root_node->setRightSib(Node::nullNode());
  root_node->setNumber(5);
  root_node->setName("root node");
  root_node->setEdgeLength(0.0);

  first_internal->setParent(root_node);
  first_internal->setLeftChild(second_internal);
  first_internal->setRightSib(Node::nullNode());
  first_internal->setNumber(4);
  first_internal->setName("first internal node");
  first_internal->setEdgeLength(0.1);

  second_internal->setParent(first_internal);
  second_internal->setLeftChild(first_leaf);
  second_internal->setRightSib(third_leaf);
  second_internal->setNumber(3);
  second_internal->setName("second internal node");
  second_internal->setEdgeLength(0.1);

  first_leaf->setParent(second_internal);
  first_leaf->setLeftChild(Node::nullNode());
  first_leaf->setRightSib(second_leaf);
  first_leaf->setNumber(0);
  first_leaf->setName("first leaf");
  first_leaf->setEdgeLength(0.1);

  second_leaf->setParent(second_internal);
  second_leaf->setLeftChild(Node::nullNode());
  second_leaf->setRightSib(Node::nullNode());
  second_leaf->setNumber(1);
  second_leaf->setName("second leaf");
  second_leaf->setEdgeLength(0.1);

  third_leaf->setParent(first_internal);
  third_leaf->setLeftChild(Node::nullNode());
  third_leaf->setRightSib(Node::nullNode());
  third_leaf->setNumber(2);
  third_leaf->setName("third leaf");
  third_leaf->setEdgeLength(0.2);

  tree->setRooted(true);
  tree->setRootNode(root_node);
  tree->setLeaves(3);

  // Note that root node is not included in preorder or levelorder
  tree->pushPreOrder(first_internal);
  tree->pushPreOrder(second_internal);
  tree->pushPreOrder(first_leaf);
  tree->pushPreOrder(second_leaf);
  tree->pushPreOrder(third_leaf);

  tree->pushLevelOrder(first_internal);
  tree->pushLevelOrder(second_internal);
  tree->pushLevelOrder(third_leaf);
  tree->pushLevelOrder(first_leaf);
  tree->pushLevelOrder(second_leaf);

  return tree;

}


std::string kpl::TreeManip::makeNewick(unsigned precision) const {

  return TreeIO::smakeNewick(getConstTree(), precision);

}




void kpl::TreeManip::refreshPreorder() {

  srefreshPreorder(getTree());

}

void kpl::TreeManip::srefreshPreorder(Tree::SharedPtr tree) {

  // Create vector of node pointers in preorder sequence

  tree->clearPreOrder();

  if (not tree->getConstRoot()) {

    return;

  }

  Node::PtrNode first_preorder = tree->getConstRoot()->getLeftChild();

  // sanity check: first preorder node should be the only child of the root node
  assert(Node::isNullNode(first_preorder->getRightSib()));

  Node::PtrNode pNode = first_preorder;

  tree->pushPreOrder(pNode);

  while (true) {

    pNode = findNextPreorder(pNode);

    if (not Node::isNullNode(pNode)) {

      tree->pushPreOrder(pNode);

    } else {

      break;

    }

  }   // end while loop

}


void kpl::TreeManip::refreshLevelorder() {

  srefreshLevelorder(getTree());

}

void kpl::TreeManip::srefreshLevelorder(Tree::SharedPtr tree) {

  if (not tree->getConstRoot()) {

    return;

  }

  // q is the buffer queue
  std::queue<Node::PtrNode > node_queue;

  // getTree()->_levelorder is the stack vector
  tree->clearLevelOrder();

  Node::PtrNode nd = tree->getConstRoot()->getLeftChild();

  // sanity check: first node should be the only child of the root node
  assert(Node::isNullNode(nd->getRightSib()));

  // Push nd onto back of queue
  node_queue.push(nd);

  while (!node_queue.empty()) {
    // pop nd off front of queue
    nd = node_queue.front();
    node_queue.pop();

    // and push it onto the stack
    tree->pushLevelOrder(nd);

    // add all children of nd to back of queue
    Node::PtrNode child = nd->getLeftChild();

    if (child) {

      node_queue.push(child);
      child = child->getRightSib();

      while (child) {

        node_queue.push(child);
        child = child->getRightSib();

      }

    }

  }   // end while loop

}


void kpl::TreeManip::renumberInternals() {

  srenumberInternals(getTree());

}

void kpl::TreeManip::srenumberInternals(Tree::SharedPtr tree) {

  assert(tree->getConstPreOrder().size() > 0);

  // Renumber internal nodes in postorder sequence
  unsigned  curr_internal = tree->numLeaves();

  for (auto nd : boost::adaptors::reverse(tree->getConstPreOrder())) {

    if (not Node::isNullNode(nd->getLeftChild())) {

      // nd is an internal node
      nd->setNumber(curr_internal);
      curr_internal++;

    }

  }

  // Root node is not included in getTree()->_preorder, so if the root node
  // is an internal node we need to number it here
  if (tree->isRooted()) {

    tree->getRootNode()->setNumber(curr_internal);
    curr_internal++;

  }

  tree->setInternals(curr_internal - tree->numLeaves());

  // If the tree has polytomies, then there are Node objects stored in
  // the getTree()->_nodes vector that have not yet been numbered. These can
  // be identified because their _number is invalid.

  assert(tree->getUnUsed().empty());

  for (; curr_internal < tree->getConstNodes().size(); curr_internal++) {

    Node::PtrNode  curr = tree->getNode(curr_internal);
    tree->pushUnused(curr);

    assert(not curr->checkValidNumber());

    curr->setNumber(curr_internal);

  }

}


void kpl::TreeManip::rerootAtNodeNumber(int node_number) {

  srerootAtNodeNumber(getTree(), node_number);

}


void kpl::TreeManip::srerootAtNodeNumber(Tree::SharedPtr tree, int node_number) {

  // Locate node having _number equal to node_number
  Node::PtrNode nd = Node::nullNode();

  for (auto curr : tree->getNodes()) {

    if (curr->getNumber() == node_number) {

      nd = curr;

      break;

    }

  }

  if (!nd) {

    throw XStrom(boost::str(boost::format("no node found with number equal to %d") % node_number));

  }

  if (nd != tree->getConstRoot()) {

    if (not Node::isNullNode(nd->getLeftChild())) {

      throw XStrom(
          boost::str(boost::format("cannot currently root trees at internal nodes (e.g. node %d)") % nd->getNumber()));

    }

    srerootAtNode(tree, nd);

  }

}


void kpl::TreeManip::rerootAtNode(Node::PtrNode prospective_root) {

  srerootAtNode(getTree(), prospective_root);

}

void kpl::TreeManip::srerootAtNode(Tree::SharedPtr tree, Node::PtrNode prospective_root) {

  Node::PtrNode a = prospective_root;
  Node::PtrNode b = prospective_root->getParent();
  Node::PtrNode c = Node::nullNode();
  Node::PtrNode d = Node::nullNode();
  Node::PtrNode p = Node::nullNode();
  a->setParent(Node::nullNode());
  double tmp_edgelen = 0.0;
  double prev_edgelen = a->getEdgeLength();

  while (not Node::isNullNode(b)) {
    // Prune node a from b
    if (a == b->getLeftChild()) {

      if (not Node::isNullNode(a->getRightSib())) {

        b->setLeftChild(a->getRightSib());
        a->setRightSib(Node::nullNode());

      } else {

        b->setLeftChild(Node::nullNode());

      }

    } else {

      c = b->getLeftChild();

      while (c->getRightSib() != a) {

        c = c->getRightSib();

      }

      d = a->getRightSib();
      c->setRightSib(d);

    }

    // Graft node b onto node a (but don't unhook node b from its parent just yet)
    if (not Node::isNullNode(a->getRightSib())) {

      c = a->getLeftChild();

      while (not Node::isNullNode(c->getRightSib())) {

        c = c->getRightSib();

      }

      c->setRightSib(b);

    } else {

      a->setLeftChild(b);

    }

    // Rotate
    p = a;
    a = b;
    b = b->getParent();
    a->setParent(p);

    // Swap nd's edge length with its new parent's edge length
    tmp_edgelen = a->getEdgeLength();
    a->setEdgeLength(prev_edgelen);
    prev_edgelen = tmp_edgelen;
  }
  prospective_root->setEdgeLength(0.0);
  tree->setRootNode(prospective_root);
  srefreshPreorder(tree);
  srefreshLevelorder(tree);

}


void kpl::TreeManip::buildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies) {

  setTree(TreeIO::buildFromNewick(newick, rooted, allow_polytomies));

}


void kpl::TreeManip::storeSplits(std::set<Split> &splitset) {

  sstoreSplits(getTree(), splitset);

}

void kpl::TreeManip::sstoreSplits(Tree::SharedPtr tree, std::set<Split> &splitset) {
  // Start by clearing and resizing all splits
  for (auto nd : tree->getNodes()) {

    nd->getSplit().resize(tree->numLeaves());

  }

  // Now do a postorder traversal and add the bit corresponding
  // to the current node in its parent node's split
  for (auto nd : boost::adaptors::reverse(tree->getConstPreOrder())) {
    if (nd->getLeftChild()) {
      // add this internal node's split to splitset
      splitset.insert(nd->getSplit());
    } else {
      // set bit corresponding to this leaf node's number
      nd->getSplit().setBitAt(nd->getNumber());
    }

    if (nd->getParent()) {
      // parent's bits are the union of the bits set in all its children
      nd->getParent()->getSplit().addSplit(nd->getSplit());

    }

  }

}


void kpl::TreeManip::selectAllPartials() {

  sselectAllPartials(getTree());

}

void kpl::TreeManip::sselectAllPartials(Tree::SharedPtr tree) {

  for (auto nd : tree->getNodes()) {

    nd->selectPartial();

  }

}


void kpl::TreeManip::deselectAllPartials() {

  sdeselectAllPartials(getTree());

}

void kpl::TreeManip::sdeselectAllPartials(Tree::SharedPtr tree) {

  for (auto nd : tree->getNodes()) {

    nd->deselectPartial();

  }

}

void kpl::TreeManip::selectAllTMatrices() {

  sselectAllTMatrices(getTree());

}


void kpl::TreeManip::sselectAllTMatrices(Tree::SharedPtr tree) {

  for (auto nd : tree->getNodes()) {

    nd->selectTMatrix();

  }

}


void kpl::TreeManip::deselectAllTMatrices() {

  sdeselectAllTMatrices(getTree());

}

void kpl::TreeManip::sdeselectAllTMatrices(Tree::SharedPtr tree) {

  for (auto nd : tree->getNodes()) {

    nd->deselectTMatrix();

  }

}


void kpl::TreeManip::selectPartialsHereToRoot(Node::PtrNode  a) {

  a->selectPartial();

  while (a->getParent()) {

    a = a->getParent();
    a->selectPartial();

  }

}


void kpl::TreeManip::flipPartialsAndTMatrices() {

  sflipPartialsAndTMatrices(getTree());

}

void kpl::TreeManip::sflipPartialsAndTMatrices(Tree::SharedPtr tree) {

  for (auto nd : tree->getNodes()) {

    if (nd->isSelPartial()) {

      if (nd->isAltPartial()) {

        nd->clearAltPartial();

      }
      else {

        nd->setAltPartial();

      }

    }

    if (nd->isSelTMatrix()) {

      if (nd->isAltTMatrix()) {

        nd->clearAltTMatrix();

      }
      else {

        nd->setAltTMatrix();

      }

    }

  }

}

void kpl::TreeManip::LargetSimonSwap(Node::PtrNode  a, Node::PtrNode  b) {

  sLargetSimonSwap(getTree(), a, b);

}


void kpl::TreeManip::sLargetSimonSwap(Tree::SharedPtr tree, Node::PtrNode  a, Node::PtrNode  b) {
  // a and b are the ends of the selected 3-edge path in a Larget-Simon move
  // The 3-edge path is indicated by parentheses around the nodes involved.
  // x is always the parent of a
  // y can be the parent of b (case 1) or the child of b (case 2)

  Node::PtrNode  x = a->getParent();
  assert(x);

  Node::PtrNode  y = x->getParent();
  assert(y);

  if (y == b->getParent()) {
    // Case 1: y is the parent of b
    //
    //    (a) d  e             (b) d  e
    //      \ | /                \ | /
    //       \|/                  \|/
    //       (x) f (b)            (x) f (a)    Swap a and b, leaving everything
    //         \ | /                \ | /      else as is
    //          \|/     ==>          \|/
    //          (y)                  (y)
    //           |                    |
    //           |                    |
    //           c                    c
    //

    // Detach a from tree
    if (a == x->getLeftChild()) {

      x->setLeftChild(a->getRightSib());

    } else {

      Node::PtrNode  child = x->getLeftChild();

      while (child->getRightSib() != a) {

        child = child->getRightSib();

      }

      child->setRightSib(a->getRightSib());

    }

    // Detach a
    a->setParent(Node::nullNode());
    a->setRightSib(Node::nullNode());

    // Detach b from tree
    if (b == y->getLeftChild()) {

      y->setLeftChild(b->getRightSib());

    } else {

      Node::PtrNode  child = y->getLeftChild();

      while (child->getRightSib() != b) {

        child = child->getRightSib();

      }
      child->setRightSib(b->getRightSib());

    }

    // Detach b
    b->setParent(Node::nullNode());
    b->setRightSib(Node::nullNode());

    // Reattach a to y as left child
    a->setRightSib(y->getLeftChild());
    y->setLeftChild(a);
    a->setParent(y);

    // Reattach b to x as left child
    b->setRightSib(x->getLeftChild());
    x->setLeftChild(b);
    b->setParent(x);

  }
  else {
    // Case 2: y is the child of b
    //
    //    (a) d  e             (a) f  c
    //      \ | /                \ | /
    //       \|/                  \|/
    //       (x) f  c            (x) d  e    swap x's children (except a)
    //         \ | /               \ | /     with y's children (except x)
    //          \|/     ==>         \|/
    //          (y)                 (y)
    //           |                   |
    //           |                   |
    //          (b)                 (b)
    assert(b == y->getParent());

    // Create a stack of x children excluding a
    std::stack<Node::PtrNode > xchildren;
    Node::PtrNode  xchild = x->getLeftChild();
    while(not Node::isNullNode(xchild)) {

      if (xchild != a) {

        xchildren.push(xchild);

      }

      xchild = xchild->getRightSib();

    }

    // Create a stack of y children excluding x
    std::stack<Node::PtrNode > ychildren;
    Node::PtrNode  ychild = y->getLeftChild();
    while(not Node::isNullNode(ychild)) {

      if (ychild != x) {

        ychildren.push(ychild);

      }

      ychild = ychild->getRightSib();

    }

  // Reattach the y children to the right of a
  x->setLeftChild(a);
  a->setParent(x);
  a->setRightSib(Node::nullNode());
  Node::PtrNode left_child = a;
  while(not ychildren.empty()) {

    Node::PtrNode child = ychildren.top();
    ychildren.pop();
    child->setParent(x);
    child->setRightSib(Node::nullNode());
    left_child->setRightSib(child);
    left_child = child;

  }

    // Reattach the x children as siblings to the right of x
    x->getParent()->setLeftChild(x);
    x->setRightSib(Node::nullNode());
    left_child = x;
    while(not xchildren.empty()) {

      Node::PtrNode child = xchildren.top();
      xchildren.pop();
      child->setParent(x->getParent());
      child->setRightSib(Node::nullNode());
      left_child->setRightSib(child);
      left_child = child;

    }


  }

  srefreshPreorder(tree);
  srefreshLevelorder(tree);

}


kpl::Node::PtrNode  kpl::TreeManip::randomInternalEdge(double uniform_deviate) {   ///begin_randomInternalEdge

  return srandomInternalEdge(getTree(), uniform_deviate);

}


/*
kpl::Node::PtrNode  kpl::TreeManip::srandomInternalEdge(Tree::SharedPtr tree, double uniform_deviate) {   ///begin_randomInternalEdge

  assert(uniform_deviate >= 0.0);
  assert(uniform_deviate < 1.0);

  // Unrooted case:                        Rooted case:
  //
  // 2     3     4     5                   1     2     3     4
  //  \   /     /     /                     \   /     /     /
  //   \ /     /     /                       \ /     /     /
  //    8     /     /                         7     /     /
  //     \   /     /                           \   /     /
  //      \ /     /                             \ /     /
  //       7     /                               6     /
  //        \   /                                 \   /
  //         \ /                                   \ /
  //          6   nleaves = 5                       5   nleaves = 4
  //          |   num_internal_edges = 2            |   num_internal_edges = 2
  //          |   choose node 7 or node 8           |   choose node 6 or node 7
  //          1                                    root
  //
  // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
  //
  // Note: _preorder is actually a vector of T *, but is shown here as a
  // vector of integers solely to illustrate the algorithm below

  int num_internal_edges = tree->getConstPreOrder().size() - tree->numLeaves() - (tree->isRooted() ? 0 : 1);
  if (num_internal_edges < 0) {
    // Star tree: return hub node, which is the first node in the preorder sequence
    return tree->getPreOrder()[0];
  }

  // Add one to skip first node in _preorder vector, which is an internal node whose edge
  // is either a terminal edge (if tree is unrooted) or invalid (if tree is rooted)
  unsigned index_of_chosen = 1 + (unsigned)std::floor(uniform_deviate*num_internal_edges);

  unsigned internal_nodes_visited = 0;
  Node::PtrNode  chosen_node = Node::nullNode();
  for (auto nd : tree->getPreOrder()) {

    if (not Node::isNullNode(nd->getLeftChild())) {

      if (internal_nodes_visited == index_of_chosen) {
        chosen_node = nd;
        break;

      }
      else {

        ++internal_nodes_visited;

      }

    }

  }
  assert(chosen_node);
  return chosen_node;
}   ///end_randomInternalEdge

*/

kpl::Node::PtrNode  kpl::TreeManip::srandomInternalEdge(Tree::SharedPtr tree, double uniform_deviate) {   ///begin_randomInternalEdge

  assert(uniform_deviate >= 0.0);
  assert(uniform_deviate < 1.0);

  // Unrooted case:                        Rooted case:
  //
  // 2     3     4     5                   1     2     3     4
  //  \   /     /     /                     \   /     /     /
  //   \ /     /     /                       \ /     /     /
  //    8     /     /                         7     /     /
  //     \   /     /                           \   /     /
  //      \ /     /                             \ /     /
  //       7     /                               6     /
  //        \   /                                 \   /
  //         \ /                                   \ /
  //          6   nleaves = 5                       5   nleaves = 4
  //          |   num_internal_edges = 2            |   num_internal_edges = 2
  //          |   choose node 7 or node 8           |   choose node 6 or node 7
  //          1                                    root
  //
  // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
  //
  // Note: _preorder is actually a vector of T *, but is shown here as a
  // vector of integers solely to illustrate the algorithm below

  int num_internal_edges = (unsigned) tree->_preorder.size() - tree->_nleaves - (tree->isRooted() ? 0 : 1);
  if (num_internal_edges < 0) {
    // Star tree: return hub node, which is the first node in the preorder sequence
    return tree->_preorder[0];
  }

  // Add one to skip first node in _preorder vector, which is an internal node whose edge
  // is either a terminal edge (if tree is unrooted) or invalid (if tree is rooted)
  unsigned index_of_chosen = 1 + (unsigned)std::floor(uniform_deviate*num_internal_edges);

  unsigned internal_nodes_visited = 0;
  Node::PtrNode  chosen_node = 0;
  for (auto nd : tree->_preorder) {
    if (nd->getLeftChild()) {
      if (internal_nodes_visited == index_of_chosen) {
        chosen_node = nd;
        break;
      }
      else
        ++internal_nodes_visited;
    }
  }
  assert(chosen_node);
  return chosen_node;
}   ///end_randomInternalEdge




/*

kpl::Node::PtrNode  kpl::TreeManip::randomInternalEdge(double uniform_deviate) {

  assert(uniform_deviate >= 0.0);
  assert(uniform_deviate < 1.0);

  // Unrooted case:                        Rooted case:
  //
  // 2     3     4     5                   1     2     3     4
  //  \   /     /     /                     \   /     /     /
  //   \ /     /     /                       \ /     /     /
  //    8     /     /                         7     /     /
  //     \   /     /                           \   /     /
  //      \ /     /                             \ /     /
  //       7     /                               6     /
  //        \   /                                 \   /
  //         \ /                                   \ /
  //          6   nleaves = 5                       5   nleaves = 4
  //          |   num_internal_edges = 2            |   num_internal_edges = 2
  //          |   choose node 7 or node 8           |   choose node 6 or node 7
  //          1                                    root
  //
  // _preorder = [6, 7, 8, 2, 3, 4, 5]     _preorder = [5, 6, 7, 1, 2, 3, 4]
  //
  // Note: _preorder is actually a vector of T *, but is shown here as a
  // vector of integers solely to illustrate the algorithm below

  unsigned num_internal_edges = (unsigned)getTree()->_preorder.size() - getTree()->_nleaves - (getTree()->_is_rooted ? 0 : 1);

  // Add one to skip first node in _preorder vector, which is an internal node whose edge
  // is either a terminal edge (if tree is unrooted) or invalid (if tree is rooted)
  unsigned index_of_chosen = 1 + (unsigned) std::floor(uniform_deviate*num_internal_edges);

  unsigned internal_nodes_visited = 0;
  Node::PtrNode  chosen_node = 0;
  for (auto nd : getTree()->_preorder) {

    if (nd->_left_child) {

      if (internal_nodes_visited == index_of_chosen) {

        chosen_node = nd;
        break;

      }
      else {

        ++internal_nodes_visited;

      }

    }

  }

  assert(chosen_node);
  return chosen_node;

}

*/

unsigned kpl::TreeManip::calcResolutionClass() const {

  return scalcResolutionClass(getConstTree());

}

unsigned kpl::TreeManip::scalcResolutionClass(Tree::ConstSharedPtr tree) {

  return tree->numInternals();

}

unsigned kpl::TreeManip::countEdges() const {

  return scountEdges(getConstTree());

}

unsigned kpl::TreeManip::scountEdges(Tree::ConstSharedPtr tree)  {

  return tree->getConstPreOrder().size();

}


unsigned kpl::TreeManip::countChildren(Node::ConstPtrNode  nd) {

  assert(nd);

  unsigned nchildren = 0;

  Node::PtrNode  child = nd->getLeftChild();

  while (child) {

    nchildren++;
    child = child->getRightSib();

  }

  return nchildren;

}


unsigned kpl::TreeManip::countInternals() const {

  return scountInternals(getConstTree()) ;

}

unsigned kpl::TreeManip::scountInternals(Tree::ConstSharedPtr tree) {

  unsigned internals = 0;

  for (auto nd : tree->getConstPreOrder()) {

    if (not Node::isNullNode(nd->getLeftChild())) {

      internals++;

    }

  }

  return internals;

}


kpl::Node::PtrNode  kpl::TreeManip::findNextPreorder(Node::ConstPtrNode  node) {

  assert(node);

  Node::PtrNode  next = Node::nullNode();

  if (Node::isNullNode(node->getLeftChild()) && Node::isNullNode(node->getRightSib())) {
    // nd has no children and no siblings, so next preorder is the right sibling of
    // the first ancestral node that has a right sibling.
    Node::PtrNode  ancestor = node->getParent();
    while (not Node::isNullNode(ancestor) && Node::isNullNode(ancestor->getRightSib())) {

      ancestor = ancestor->getParent();

    }
    if (not Node::isNullNode(ancestor)) {
      // We found an ancestor with a right sibling
      next = ancestor->getRightSib();
    }
    else {
      // nd is last preorder node in the tree
      next = Node::nullNode();
    }
  }
  else if (not Node::isNullNode(node->getRightSib()) && Node::isNullNode(node->getLeftChild())) {
    // nd has no children (it is a tip), but does have a sibling on its right
    next = node->getRightSib();
  }
  else if (not Node::isNullNode(node->getLeftChild()) && Node::isNullNode(node->getRightSib())) {
    // nd has children (it is an internal node) but no siblings on its right
    next = node->getLeftChild();
  }
  else {
    // nd has both children and siblings on its right
    next = node->getLeftChild();
  }

  return next;

}



kpl::Node::PtrNode  kpl::TreeManip::findLeftSib(Node::ConstPtrNode  node) {

  assert(node);
  assert(node->getParent());

  Node::PtrNode  child = node->getParent()->getLeftChild();

  while (not Node::isNullNode(child) && child->getRightSib() != node) {

    child = child->getRightSib();

  }

  return child;

}


kpl::Node::PtrNode  kpl::TreeManip::findRightmostChild(Node::ConstPtrNode  pNode) {

  assert(pNode);

  Node::PtrNode  child = pNode->getLeftChild();

  while (child->getRightSib()) {

    child = child->getRightSib();

  }

  return child;

}


kpl::Node::PtrNode  kpl::TreeManip::findLastPreorderInClade(Node::PtrNode  start) {

  assert(start);

  Node::PtrNode  curr = start;
  Node::PtrNode  rchild = findRightmostChild(curr);

  while (not Node::isNullNode(rchild)) {

    curr = rchild;
    rchild = findRightmostChild(curr);

  }

  return curr;

}


void kpl::TreeManip::insertSubtreeOnLeft(Node::PtrNode  subtree, Node::PtrNode  parent) {

  assert(parent);
  assert(subtree);
  subtree->setRightSib(parent->getLeftChild());
  subtree->setParent(parent);
  parent->setLeftChild(subtree);

}


void kpl::TreeManip::insertSubtreeOnRight(Node::PtrNode  subtree, Node::PtrNode  parent) {

  assert(parent);
  assert(subtree);

  subtree->setRightSib(Node::nullNode());
  subtree->setParent(parent);

  if (not Node::isNullNode(parent->getLeftChild())) {

    Node::PtrNode  u_rchild = findRightmostChild(parent);
    u_rchild->setRightSib(subtree);

  }
  else {

    parent->setLeftChild(subtree);

  }

}


void kpl::TreeManip::detachSubtree(Node::PtrNode  subtree) {

  assert(subtree);
  assert(not Node::isNullNode(subtree->getParent()));

  // Save pointers to relevant nodes
  Node::PtrNode  subtree_leftsib  = findLeftSib(subtree);
  Node::PtrNode  subtree_rightsib = subtree->getRightSib();
  Node::PtrNode  subtree_parent   = subtree->getParent();

  // Completely detach s and seal up the wound
  subtree->setParent(Node::nullNode());
  subtree->setRightSib(Node::nullNode());

  if (not Node::isNullNode(subtree_leftsib)) {

    subtree_leftsib->setRightSib(subtree_rightsib);

  }
  else {

    subtree_parent->setLeftChild(subtree_rightsib);

  }

}

void kpl::TreeManip::rectifyNumInternals(int incr) {

  srectifyNumInternals(getTree(), incr);

}

void kpl::TreeManip::srectifyNumInternals(Tree::SharedPtr tree, int incr) {

  assert(tree->getConstNodes().size() == tree->getUnUsed().size() + tree->numLeaves() + tree->numInternals() + incr);
  tree->setInternals(incr + tree->numInternals());

}


void kpl::TreeManip::refreshNavigationPointers() {

  srefreshNavigationPointers(getTree());

}

void kpl::TreeManip::srefreshNavigationPointers(Tree::SharedPtr tree) {

  srefreshPreorder(tree);
  srefreshLevelorder(tree);

}


kpl::Node::PtrNode  kpl::TreeManip::getUnusedNode() {

  return sgetUnusedNode(getTree());

}

kpl::Node::PtrNode  kpl::TreeManip::sgetUnusedNode(Tree::SharedPtr tree) {

  assert(tree->getUnUsed().size() > 0);

  Node::PtrNode  node = tree->getUnUsed().top();
  tree ->popUnused();
  node->clearPointers();

  return node;

}

void kpl::TreeManip::putUnusedNode(Node::PtrNode  node) {

  sputUnusedNode(getTree(), node);

}

void kpl::TreeManip::sputUnusedNode(Tree::SharedPtr tree, Node::PtrNode  node) {

  tree->pushUnused(node);

}