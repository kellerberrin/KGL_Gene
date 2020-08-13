//
// Created by kellerberrin on 18/12/19.
//

#include "kpl_mcmc_treeupdater.h"
#include "kpl_strom.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::TreeUpdater::TreeUpdater() {
  // std::cout << "Creating a TreeUpdater" << std::endl;
  Updater::clear();
  name("Tree Topology and Edge Proportions");
  reset();
}


kpl::TreeUpdater::~TreeUpdater() {
  // std::cout << "Destroying a TreeUpdater" << std::endl;
}


void kpl::TreeUpdater::reset() {

  _topology_changed       = false;
  _orig_edgelen_top       = 0.0;
  _orig_edgelen_middle    = 0.0;
  _orig_edgelen_bottom    = 0.0;
  logHastingsRatio(0.0);
  _case                   = 0;
  _x                      = nullptr;
  _y                      = nullptr;
  _a                      = nullptr;
  _b                      = nullptr;

}


void kpl::TreeUpdater::starTreeMove() {
  // Choose focal 2-edge segment to modify
  _orig_edgelen_middle = 0.0;

  // Choose the first edge
  _a = chooseRandomChild(_x, 0, false);
  _orig_edgelen_top = _a->getEdgeLength();

  // Choose the second edge
  _b = chooseRandomChild(_x, _a, true);

  if (!_b) {

    _b = _x;

  }
  _orig_edgelen_bottom = _b->getEdgeLength();

  // Note that _a must be a child of _x, but _b may either be a different child of _x or _x itself
  double u = lot()->uniform();
  double new_edgelen_top    = u*(_orig_edgelen_top + _orig_edgelen_bottom);
  double new_edgelen_bottom = (1.0 - u)*(_orig_edgelen_top + _orig_edgelen_bottom);

  // Hastings ratio and Jacobian are both 1 under Gamma-Dirichlet parameterization
  logHastingsRatio(0.0);
  logJacobian(0.0);

  // Change edge lengths and flag partials and transition matrices for recalculation
  treeManipulator()->selectPartialsHereToRoot(_x);
  _a->setEdgeLength(new_edgelen_top);
  _a->selectTMatrix();
  _b->setEdgeLength(new_edgelen_bottom);
  _b->selectTMatrix();

}


kpl::Node::PtrNode  kpl::TreeUpdater::chooseRandomChild(Node::PtrNode  x, Node::PtrNode  avoid, bool parent_included) {
  // Count number of children of x
  unsigned n = 0;
  Node::PtrNode  child = x->getLeftChild();

  while (child) {

    if (child != avoid) {

      n++;

    }
    child = child->getRightSib();

  }

  // Choose random child index
  unsigned upper = n + (parent_included ? 1 : 0);
  unsigned chosen = lot()->randint(0,upper - 1);

  // If chosen < n, then find and return that particular child
  if (chosen < n) {

    child = x->getLeftChild();
    unsigned i = 0;

    while (child) {

      if (child != avoid && i == chosen) {

        return child;

      }
      else if (child != avoid) {

        i++;

      }

      child = child->getRightSib();

    }

  }

  // If chosen equals n, then the parent was chosen, indicated by returning NULL
  return NULL;

}


double kpl::TreeUpdater::calcLogPrior() {

  double log_topology_prior = calcLogTopologyPrior();
  double log_edge_length_prior = Updater::calcEdgeLengthPrior();
  return log_topology_prior + log_edge_length_prior;

}


void kpl::TreeUpdater::proposeNewState() {

  _case = 0;
  _topology_changed = false;

//  ExecEnv::log().info("TreeUpdater::proposeNewState() {}", treeManipulator()->getTree()->treeDescription());
//  if (treeManipulator()->getTree()->isRooted()) {

//    ExecEnv::log().info("TreeUpdater::proposeNewState() Tree is unrooted");
//    treeManipulator()->getTree()->setRoot(Node::nullNode());
//    ExecEnv::log().info("TreeUpdater::proposeNewState() {}", treeManipulator()->getTree()->treeDescription());

//  }

  assert(!treeManipulator()->getTree()->isRooted());

  // Choose random internal node x that is not the root and has parent y that is also not the root.
  // After choosing x (and thus y), choose a and b randomly to create a contiguous 3-edge segment.
  //
  //        a
  //   \ | /
  //    \|/
  //     x
  //      \ | /
  //       \|/
  //        y
  //        |
  //        |
  //        b
  // For the star tree, there is only one internal node. In this case, only choose
  // two edges and modify them (no change in tree topology is possible)
  //
  //           a
  //      \ | /
  //       \|/
  //        x
  //        |
  //        |
  //        b
  //

  _x = treeManipulator()->randomInternalEdge(lot()->uniform());
  _orig_edgelen_middle = _x->getEdgeLength();

  // The only child of the root node will be chosen only if the tree equals the star tree
  // in which case we want to perform a starTreeMove rather than Larget-Simon
  _star_tree_move = false;
  if (_x->getParent() && !_x->getParent()->getParent()) {
    _star_tree_move = true;
    starTreeMove();
    return;
  }

  _y = _x->getParent();

  // Choose focal 3-edge segment to modify
  // Begin by randomly choosing one child of x to be node _a
  _a = chooseRandomChild(_x, 0, false);
  _orig_edgelen_top = _a->getEdgeLength();

  // Now choose a different child of x (or the parent) to be node _b
  _b = chooseRandomChild(_y, _x, true);
  bool b_is_child_of_y;
  if (_b) {

    b_is_child_of_y = true;
    _orig_edgelen_bottom = _b->getEdgeLength();

  }
  else {

    b_is_child_of_y = false;
    _b = _y->getParent();
    _orig_edgelen_bottom = _y->getEdgeLength();

  }

  // Symmetric move: Hastings ratio = 1
  logHastingsRatio(0.0);

  // Decide where along focal path (starting from top) to place moved node
  double new_focal_path_length = _orig_edgelen_top + _orig_edgelen_middle + _orig_edgelen_bottom;
  double u = lot()->uniform();
  double new_attachment_point = u*new_focal_path_length;

  if (new_attachment_point <= Node::smallestEdgeLength()) {

    new_attachment_point = Node::smallestEdgeLength();

  }
  else if (new_focal_path_length - new_attachment_point <= Node::smallestEdgeLength()) {

    new_attachment_point = new_focal_path_length - Node::smallestEdgeLength();

  }

  // Decide which node(s) to move, and whether the move involves a topology change
  u = lot()->uniform();
  bool x_node_slides = (bool)(u < 0.5);

  if (x_node_slides) {

    // _x slides toward _y
    _topology_changed = (new_attachment_point > _orig_edgelen_top + _orig_edgelen_middle);
    if (_topology_changed) {

      treeManipulator()->LargetSimonSwap(_a, _b);
      if (b_is_child_of_y) {

        // LargetSimonSwap case 1: a swapped with b
        _a->setEdgeLength(new_focal_path_length - new_attachment_point);
        _x->setEdgeLength(new_attachment_point - _orig_edgelen_top - _orig_edgelen_middle);
        _b->setEdgeLength(_orig_edgelen_top + _orig_edgelen_middle);
        _case = 1;

      } else {

        // LargetSimonSwap case 2: x's children (except a) swapped with y's children (except b)
        _a->setEdgeLength(_orig_edgelen_top + _orig_edgelen_middle);
        _x->setEdgeLength(new_attachment_point - _orig_edgelen_top - _orig_edgelen_middle);
        _y->setEdgeLength(new_focal_path_length - new_attachment_point);
        _case = 2;

      }

    } else {

      _a->setEdgeLength(new_attachment_point);
      _x->setEdgeLength(_orig_edgelen_top + _orig_edgelen_middle - new_attachment_point);

      if (b_is_child_of_y) {

        _b->setEdgeLength(_orig_edgelen_bottom);
        _case = 3;

      } else {

        _y->setEdgeLength(_orig_edgelen_bottom);
        _case = 4;

      }

    }

  } else {

    // _y slides toward _x
    _topology_changed = (new_attachment_point < _orig_edgelen_top);
    if (_topology_changed) {

      treeManipulator()->LargetSimonSwap(_a, _b);

      if (b_is_child_of_y) {

        // LargetSimonSwap case 1: a swapped with b
        _a->setEdgeLength(_orig_edgelen_middle + _orig_edgelen_bottom);
        _x->setEdgeLength(_orig_edgelen_top - new_attachment_point);
        _b->setEdgeLength(new_attachment_point);
        _case = 5;

      } else {

        // LargetSimonSwap case 2: x's children (except a) swapped with y's children (except b)
        _a->setEdgeLength(new_attachment_point);
        _x->setEdgeLength(_orig_edgelen_top - new_attachment_point);
        _y->setEdgeLength(_orig_edgelen_middle + _orig_edgelen_bottom);
        _case = 6;

      }

    } else {

      _a->setEdgeLength(_orig_edgelen_top);
      _x->setEdgeLength(new_attachment_point - _orig_edgelen_top);

      if (b_is_child_of_y) {

        _b->setEdgeLength(new_focal_path_length - new_attachment_point);
        _case = 7;

      } else {

        _y->setEdgeLength(new_focal_path_length - new_attachment_point);
        _case = 8;

      }

    }

  }

  // flag partials and transition matrices for recalculation
  treeManipulator()->selectPartialsHereToRoot(_x);
  _a->selectTMatrix();
  _x->selectTMatrix();

  if (_case == 2 || _case == 4 || _case == 6 || _case == 8) {
    // In these cases b is below y, so it is y's edge that is modified
    _y->selectTMatrix();
  } else {
    // In these cases b is above y, so it is b's edge that is modified
    _b->selectTMatrix();
  }

}


void kpl::TreeUpdater::revert() {

  if (_star_tree_move) {

    _a->setEdgeLength(_orig_edgelen_top);
    _b->setEdgeLength(_orig_edgelen_bottom);

  }
  else {

    assert(_case > 0 && _case < 9);

    if (_case == 2 || _case == 6) {

      treeManipulator()->LargetSimonSwap(_a, _b);

    }
    else if (_case == 1 || _case == 5) {

      treeManipulator()->LargetSimonSwap(_b, _a);

    }

    _a->setEdgeLength(_orig_edgelen_top);
    _x->setEdgeLength(_orig_edgelen_middle);

    if (_case == 1 || _case == 3 || _case == 5 || _case == 7) {

      _b->setEdgeLength(_orig_edgelen_bottom);

    }
    else {

      _y->setEdgeLength(_orig_edgelen_bottom);

    }

  }

}

