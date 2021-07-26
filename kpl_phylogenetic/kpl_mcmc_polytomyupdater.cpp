//
// Created by kellerberrin on 20/12/19.
//

#include "kpl_mcmc_polytomyupdater.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::PolytomyUpdater::PolytomyUpdater() {
  // std::cout << "Creating a PolytomyUpdater" << std::endl;
  Updater::clear();
  name("Polytomies");
  reset();

}


kpl::PolytomyUpdater::~PolytomyUpdater() {
  // std::cout << "Destroying a PolytomyUpdater" << std::endl;

}


void kpl::PolytomyUpdater::reset() {

  _tree_length          = 0.0;
  _new_edge_proportion  = 0.0;
  _orig_edge_proportion = 0.0;
  _phi                  = 0.5;
  _orig_par             = 0;
  _orig_lchild          = 0;
  _polytomy_size        = 0;
  _num_polytomies       = 0;
  _add_edge_proposed    = false;
  _polytomies.clear();

}

double kpl::PolytomyUpdater::calcLogPrior() {

  double log_prior = 0.0;
  log_prior += Updater::calcLogTopologyPrior();
// Note that this function was named "double calcLogEdgeLengthPrior() const" in the source code that
// was with the tutorial. Some source code in kpl_mcmc_polytomyupdater.cc referred to this alternate
// function name in the derived class function "double kpl::PolytomyUpdater::calcLogPrior()" and was modified.
// to use "double calcEdgeLengthPrior() const" instead.
  log_prior += Updater::calcEdgeLengthPrior();

  return log_prior;

}


void kpl::PolytomyUpdater::proposeNewState() {

  std::shared_ptr<Tree> tree = treeManipulator()->getTree();

  // Translate tuning parameter lambda_ into the maximum possible proportion
  // that a newly created edge could have
  if (lambda() > 1000.0) {

    lambda(1000.0);

  }
  else if (lambda() < 1.0) {

    lambda(1.0);

  }

  _phi = 1.0 - std::exp(- lambda());

  // Compute number of internal nodes in a fully resolved tree
  unsigned num_internals_in_fully_resolved_tree = 0;

  if (tree->isRooted()) {

    num_internals_in_fully_resolved_tree = tree->numLeaves();

  }
  else {

    num_internals_in_fully_resolved_tree = tree->numLeaves() - 2;

  }

  // Compute tree length before proposed move
  _tree_length = treeManipulator()->calcTreeLength();

  // Determine whether starting tree is fully resolved or the star tree
  unsigned num_internals_before = tree->numInternals();
  unsigned num_leaves_before = tree->numLeaves();
  unsigned num_internal_edges_before = num_internals_before - (tree->isRooted() ? 2 : 1);
  unsigned total_edges_before = num_leaves_before + num_internal_edges_before;
  bool fully_resolved_before = (num_internals_in_fully_resolved_tree == num_internals_before);
  bool star_tree_before = (tree->numInternals() == 1);

  // Rebuild _polytomies vector, which is a vector of pointers to Node objects having more than 2 children
  refreshPolytomies();
  _num_polytomies = (unsigned)_polytomies.size();

  // Determine whether an add edge move is proposed (alternative is a delete edge move)
  if (star_tree_before) {

    _add_edge_proposed = true;

  }
  else if (fully_resolved_before) {

    _add_edge_proposed = false;

  }
  else {

    _add_edge_proposed = (lot()->uniform() < 0.5);

  }

  if (_add_edge_proposed) {
    // Choose a polytomy at random to split

    unsigned i = lot()->randint(0, _num_polytomies-1);

    Node::PtrNode pNode = _polytomies[i];

    _polytomy_size = 1 + treeManipulator()->countChildren(pNode);

    // Add an edge to split up polytomy at nd, moving a random subset
    // of the spokes to the (new) left child of nd
    proposeAddEdgeMove(pNode);
    Node::PtrNode  new_nd = pNode->getLeftChild();

    double TL_after_add_edge = treeManipulator()->calcTreeLength();
    assert(std::fabs(_tree_length - TL_after_add_edge) < 1.e-8);

    // Compute the log of the Hastings ratio
    double log_hastings_ratio  = 0.0;
    log_hastings_ratio += std::log(_num_polytomies);
    log_hastings_ratio += std::log(pow(2.0, _polytomy_size - 1) - _polytomy_size - 1);
    log_hastings_ratio -= std::log(num_internal_edges_before + 1);

    // Now multiply by the value of the quantity labeled gamma_b in the Lewis-Holder-Holsinger (2005) paper
    unsigned num_internals_after = tree->numInternals();
    assert(num_internals_after == num_internals_before + 1);
    const bool fully_resolved_after = (num_internals_after == num_internals_in_fully_resolved_tree);
    if (star_tree_before && !fully_resolved_after) {

      log_hastings_ratio -= log(2.0);

    }
    else if (fully_resolved_after && !star_tree_before) {

      log_hastings_ratio += log(2.0);

    }

    // Compute the log of the Jacobian
    double log_jacobian = 0.0;
    log_jacobian += std::log(_phi);
    log_jacobian += (total_edges_before - 1)*std::log(1.0 - _new_edge_proportion);

    // update the mcmc
    logHastingsRatio(log_hastings_ratio);
    logJacobian(log_jacobian);

    // flag partials and transition_ matrices for recalculation
    treeManipulator()->selectPartialsHereToRoot(new_nd);
    new_nd->selectTMatrix();
  }
  else {
    // Choose an internal edge at random and delete it to create a polytomy
    // (or a bigger polytomy if there is already a polytomy)
    double select_rand = lot()->uniform();

    Node::PtrNode pNode = treeManipulator()->randomInternalEdge(select_rand);

    proposeDeleteEdgeMove(pNode);

    double TL_after_del_edge = treeManipulator()->calcTreeLength();
    assert(std::fabs(_tree_length - TL_after_del_edge) < 1.e-8);

    // Compute the log of the Hastings ratio
    double log_hastings_ratio  = 0.0;
    log_hastings_ratio += std::log(num_internal_edges_before);
    log_hastings_ratio -= std::log(_num_polytomies);
    log_hastings_ratio -= std::log(pow(2.0, _polytomy_size - 1) - _polytomy_size - 1);

    // Now multiply by the value of the quantity labeled gamma_b in the Lewis-Holder-Holsinger (2005) paper
    // Now multiply by the value of the quantity labeled gamma_d in the paper
    unsigned num_internals_after = tree->numInternals();
    assert(num_internals_after == num_internals_before - 1);
    const bool star_tree_after = (num_internals_after == (tree->isRooted() ? 2 : 1));
    if (fully_resolved_before && !star_tree_after) {

      log_hastings_ratio -= log(2.0);

    }
    else if (star_tree_after && !fully_resolved_before) {

      log_hastings_ratio += log(2.0);

    }

    // Compute the log Jacobian
    double log_jacobian = 0.0;
    log_jacobian -= std::log(_phi);
    log_jacobian -= (total_edges_before - 2)*std::log(1.0 - _orig_edge_proportion);

    // update the mcmc
    logHastingsRatio(log_hastings_ratio);
    logJacobian(log_jacobian);

    // flag partials and transition_ matrices for recalculation
    treeManipulator()->selectPartialsHereToRoot(pNode);

  }

}



void kpl::PolytomyUpdater::proposeAddEdgeMove(Node::PtrNode  u) {

  std::shared_ptr<Tree> tree = treeManipulator()->getTree();

  // Split up the polytomy at `u' by creating a new internal node v and a new edge
  // connecting u with v. Node u is saved as _orig_par and node v is saved
  // as _orig_lchild in case we need to revert the proposed move.
  assert(u);

  // Calculate (if necessary) the probability of each possible partitioning of the chosen polytomy
  // Select number of spokes to move over to new node
  // Note that 0 and 1 are not allowed because they
  // would leave the tree in an invalid FSM_State
  const _partition_vect_t & prob_n = computePolytomyDistribution(_polytomy_size);
  double p = lot()->uniform();
  double cum = 0.0;
  unsigned k = 0;
  bool found = false;
  for (; k < prob_n.size(); ++k) {

    double prob_k_given_n = prob_n[k];
    cum += prob_k_given_n;

    if (p < cum) {

      found = true;
      break;

    }

  }
  assert(found);
  k += 2;

  // Create the new node that will receive the k randomly-chosen spokes
  Node::PtrNode  v = treeManipulator()->getUnusedNode();
  _new_edge_proportion = lot()->uniform()*_phi;
  treeManipulator()->scaleAllEdgeLengths(1.0 - _new_edge_proportion);
  v->setEdgeLength(_new_edge_proportion*_tree_length);
  treeManipulator()->rectifyNumInternals(+1);
  treeManipulator()->insertSubtreeOnLeft(v, u);
  assert(u->getLeftChild() == v);

  // Save u and v. If revert is necessary, all of orig_lchild's nodes will be returned
  // to orig_par, and orig_lchild will be deleted.
  _orig_par    = u;
  _orig_lchild = v;

  // After the move, either v should have k spokes and the
  // other node _polytomy_size - k spokes (u and v will each
  // have 1 additional connector spoke).
  //
  // Choose k spokes randomly out of the _polytomy_size available.
  // If u->par is included, let u retain the k spokes and move
  // _polytomy_size - k spokes to v. Otherwise, move the
  // k spokes to v leaving _polytomy_size - k spokes behind.

  // Create vector of valid spokes (parent and all children of u except v)
  std::vector<Node::PtrNode > uspokes;
  uspokes.push_back(u->getParent());
  for (Node::PtrNode  uchild = u->getLeftChild(); uchild; uchild = uchild->getRightSib()) {

    if (uchild != v) {

      uspokes.push_back(uchild);

    }
  }
  assert (uspokes.size() == _polytomy_size);

  bool reverse_polarity = false;
  std::vector<Node::PtrNode > vspokes;
  typedef std::vector<Node::PtrNode >::iterator::difference_type vec_it_diff;
  for (unsigned i = 0; i < k; ++i) {

    unsigned num_u_spokes = (unsigned)uspokes.size();
    assert(num_u_spokes > 0);
    unsigned j = lot()->randint(0, num_u_spokes-1);
    Node::PtrNode  s = uspokes[j];
    if (s == u->getParent()) {

      reverse_polarity = true;

    }
    vspokes.push_back(s);
    uspokes.erase(uspokes.begin() + (vec_it_diff)j);

  }
  assert(uspokes.size() + vspokes.size() == _polytomy_size);

  if (reverse_polarity) {
    // transfer nodes in uspokes to v
    for (auto s = uspokes.begin(); s != uspokes.end(); ++s) {

      treeManipulator()->detachSubtree(*s);
      treeManipulator()->insertSubtreeOnRight(*s, v);

    }

  }
  else {
    // transfer nodes in vspokes to v
    for (auto s = vspokes.begin(); s != vspokes.end(); ++s) {

      treeManipulator()->detachSubtree(*s);
      treeManipulator()->insertSubtreeOnRight(*s, v);

    }

  }

  treeManipulator()->refreshNavigationPointers();

}


void kpl::PolytomyUpdater::proposeDeleteEdgeMove(Node::PtrNode  u) {
  // Delete the edge associated with `u' to create a polytomy (or a bigger polytomy if `u->par' was already a polytomy).
  // The supplied node u should not be the only child of the root node.
  //
  //       b       c
  //        \     /
  //         \   /
  //          \ /
  //   a       u           a   b   c
  //    \     /            \  |  /
  //     \   /              \ | /
  //      \ /                \|/
  //       v                  v
  //      /                  /
  //
  //     Before           After
  //
  // Returns the number of polytomies in the tree after the proposed delete-edge move.
  // The return value will be incorrect if the polytomies vector is not up-to-date.

  // Save u's edge length in case we need to revert
  _orig_edge_proportion = u->getEdgeLength()/_tree_length;

  // This operation should not leave the root node (which is a tip) with more than one
  // child, so check to make sure that the supplied node is not the root nor a child of root
  _orig_par = u->getParent();
  assert(_orig_par);
  assert(_orig_par->getParent());

  // Compute size of polytomy after the delete-edge move, a quantity that is needed for computing the Hastings ratio.
  // Note that one of v's children (i.e. u) is deleted but this is made up for by considering v->par, which is
  // also a spoke that counts.
  unsigned u_children = treeManipulator()->countChildren(u);
  unsigned v_children = treeManipulator()->countChildren(_orig_par);
  _polytomy_size = v_children + u_children;

  bool u_polytomy_before = (u_children > 2);
  bool v_polytomy_before = (v_children > 2);
  if (u_polytomy_before && v_polytomy_before) {
    // No. polytomies will decrease by one as a result of this delete-edge move
    --_num_polytomies;

  }
  else if (!u_polytomy_before && !v_polytomy_before) {
    // No. polytomies will increase by one as a result of this delete-edge move
    ++_num_polytomies;

  }

  // Make all of u's children siblings of u (i.e. children of u->par)
  _orig_lchild = u->getLeftChild();
  while (u->getLeftChild() != NULL) {

    Node::PtrNode  tmp = u->getLeftChild();
    treeManipulator()->detachSubtree(tmp);
    treeManipulator()->insertSubtreeOnRight(tmp, _orig_par);

  }

  treeManipulator()->detachSubtree(u);
  treeManipulator()->putUnusedNode(u);
  treeManipulator()->rectifyNumInternals(-1);

  treeManipulator()->refreshNavigationPointers();
  assert(_orig_edge_proportion < 1.0);
  double scaler = 1.0/(1.0 - _orig_edge_proportion);
  treeManipulator()->scaleAllEdgeLengths(scaler);

}


kpl::PolytomyUpdater::_partition_vect_t& kpl::PolytomyUpdater::computePolytomyDistribution(unsigned nspokes) {
  assert(nspokes > 2);

  // Only compute it if it isn't already stored in the _poly_prob map
  auto iter = _poly_prob.find(nspokes);

  if (iter == _poly_prob.end()) {
    // There is no existing probability distribution vector corresponding to nspokes
    double ln_denom = std::log(pow(2.0,nspokes-1) - nspokes - 1.0);
    _partition_vect_t v(nspokes - 3);
    unsigned first = 2;
    unsigned last = nspokes/2;
    bool nspokes_even = nspokes % 2 == 0;
    double total_prob = 0.0;
    for (unsigned x = first; x <= last; ++x) {

      double ln_numer = std::lgamma(nspokes + 1) - std::lgamma(x + 1) - std::lgamma(nspokes - x + 1);

      if (nspokes_even && x == last) {

        ln_numer -= std::log(2);

      }
      double prob_x = exp(ln_numer - ln_denom);

      if (prob_x > 1.0) {

        prob_x = 1.0;

      }

      total_prob += prob_x;
      v[x-first] = prob_x;

    }

    assert(std::fabs(total_prob - 1.0) < 1.e-8);
    _poly_prob[nspokes] = v;

  }

  return _poly_prob[nspokes];

}


void kpl::PolytomyUpdater::revert() {

  if (_add_edge_proposed) {
    // Return all of _orig_lchild's child nodes to _orig_par, then return orig_lchild to storage
    Node::PtrNode  child = _orig_lchild->getLeftChild();

    while (child) {

      Node::PtrNode  rsib = child->getRightSib();
      treeManipulator()->detachSubtree(child);
      treeManipulator()->insertSubtreeOnRight(child, _orig_par);
      child = rsib;

    }

    treeManipulator()->detachSubtree(_orig_lchild);
    treeManipulator()->putUnusedNode(_orig_lchild);
    treeManipulator()->rectifyNumInternals(-1);

    treeManipulator()->refreshNavigationPointers();
    assert(_new_edge_proportion < 1.0);
    double scaler = 1.0/(1.0 - _new_edge_proportion);
    treeManipulator()->scaleAllEdgeLengths(scaler);
    double TL_after_revert_add_edge = treeManipulator()->calcTreeLength();
    assert(std::fabs(_tree_length -TL_after_revert_add_edge) < 1.e-8);

  }
  else {

    Node::PtrNode  v = treeManipulator()->getUnusedNode();
    treeManipulator()->rectifyNumInternals(+1);
    treeManipulator()->scaleAllEdgeLengths(1.0 - _orig_edge_proportion);
    v->setEdgeLength(_orig_edge_proportion*_tree_length);
    for (Node::PtrNode  child = _orig_lchild; child;) {

      Node::PtrNode  child_rsib = child->getRightSib();
      treeManipulator()->detachSubtree(child);
      treeManipulator()->insertSubtreeOnRight(child, v);
      child = child_rsib;

    }

    treeManipulator()->insertSubtreeOnRight(v, _orig_par);

    treeManipulator()->refreshNavigationPointers();
    double TL_after_revert_del_edge = treeManipulator()->calcTreeLength();
    assert(std::fabs(_tree_length - TL_after_revert_del_edge) < 1.e-8);

  }

}


void kpl::PolytomyUpdater::refreshPolytomies() {

  std::shared_ptr<Tree> tree = treeManipulator()->getTree();
  _polytomies.clear();

  for (auto nd : tree->getConstPreOrder()) {

    unsigned s = treeManipulator()->countChildren(nd);

    if (s > 2) {

      _polytomies.push_back(nd);

    }

  }

}
