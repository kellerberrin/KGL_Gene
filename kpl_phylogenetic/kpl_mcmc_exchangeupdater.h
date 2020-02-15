//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_EXCHANGEUPDATER_H
#define KPL_MCMC_EXCHANGEUPDATER_H


#include "kpl_mcmc_dirichlet.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


class ExchangeabilityUpdater : public DirichletUpdater {

public:
  typedef std::shared_ptr< ExchangeabilityUpdater > SharedPtr;

  explicit ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix);
  ~ExchangeabilityUpdater() override;

private:

  void pullFromModel() override;
  void pushToModel() override;

  QMatrix::SharedPtr              _qmatrix;

};


} // end namespace


#endif //KPL_MCMC_EXCHANGEUPDATER_H
