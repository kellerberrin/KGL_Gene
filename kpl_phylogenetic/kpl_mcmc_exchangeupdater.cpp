//
// Created by kellerberrin on 17/12/19.
//

#include "kpl_mcmc_exchangeupdater.h"

namespace kpl = kellerberrin::phylogenetic;


kpl::ExchangeabilityUpdater::ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix) {
  // std::cout << "Creating an ExchangeabilityUpdater" << std::endl;
  DirichletUpdater::clear();
  name("Exchangeabilities");
  assert(qmatrix);
  _qmatrix = qmatrix;

}


kpl::ExchangeabilityUpdater::~ExchangeabilityUpdater() {
  // std::cout << "Destroying an ExchangeabilityUpdater" << std::endl;
}


void kpl::ExchangeabilityUpdater::pullFromModel() {

  QMatrix::freq_xchg_ptr_t xchg = _qmatrix->getExchangeabilitiesSharedPtr();
  _curr_point.assign(xchg->begin(), xchg->end());

}


void kpl::ExchangeabilityUpdater::pushToModel() {

  _qmatrix->setExchangeabilities(_curr_point);

}