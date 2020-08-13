//
// Created by kellerberrin on 18/12/19.
//

#include "kpl_mcmc_treelength.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::TreeLengthUpdater::TreeLengthUpdater() {
  // std::cout << "Creating a TreeLengthUpdater..." << std::endl;
  clear();
  name("Tree Length");

}

kpl::TreeLengthUpdater::~TreeLengthUpdater() {
  // std::cout << "Destroying a TreeLengthUpdater..." << std::endl;

}


void kpl::TreeLengthUpdater::clear() {

  Updater::clear();
  _prev_point     = 0.0;
  _curr_point     = 0.0;
  reset();

}


double kpl::TreeLengthUpdater::calcLogPrior() {

  return Updater::calcEdgeLengthPrior();

}


void kpl::TreeLengthUpdater::pullFromModel() {

  _curr_point = treeManipulator()->calcTreeLength();

}


void kpl::TreeLengthUpdater::pushToModel() {

  double scaler = _curr_point/_prev_point;
  treeManipulator()->scaleAllEdgeLengths(scaler);

}


void kpl::TreeLengthUpdater::proposeNewState() {

  // Save copy of _curr_point in case revert is necessary.
  pullFromModel();
  _prev_point = _curr_point;

  // Let _curr_point be proposed value
  double m = exp(lambda() *(lot()->uniform() - 0.5));
  _curr_point = m*_prev_point;
  pushToModel();

  // calculate log of Hastings ratio under GammaDir parameterization
  logHastingsRatio(log(m));

  // This proposal invalidates all transition matrices and partials
  treeManipulator()->selectAllPartials();
  treeManipulator()->selectAllTMatrices();

}


void kpl::TreeLengthUpdater::revert() {

  // swap _curr_point and _prev_point so that edge length scaler
  // in pushCurrentStateToModel will be correctly calculated
  double tmp = _curr_point;
  _curr_point = _prev_point;
  _prev_point = tmp;
  pushToModel();

}