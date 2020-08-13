//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_statefrequpdater.h"


namespace kpl = kellerberrin::phylogenetic;


kpl::StateFreqUpdater::StateFreqUpdater(QMatrix::SharedPtr qmatrix) {
  //std::cout << "Creating a StateFreqUpdater" << std::endl;
  DirichletUpdater::clear();
  name("State Frequencies");
  assert(qmatrix);
  _qmatrix = qmatrix;

}


kpl::StateFreqUpdater::~StateFreqUpdater() {
  //std::cout << "Destroying a StateFreqUpdater" << std::endl;
}


void kpl::StateFreqUpdater::pullFromModel() {

  QMatrix::freq_xchg_ptr_t freqs = _qmatrix->getStateFreqsSharedPtr();
  _curr_point.assign(freqs->begin(), freqs->end());

}


void kpl::StateFreqUpdater::pushToModel() {

  _qmatrix->setStateFreqs(_curr_point);

}