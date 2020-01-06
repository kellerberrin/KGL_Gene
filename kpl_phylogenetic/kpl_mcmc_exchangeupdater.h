//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_EXCHANGEUPDATER_H
#define KPL_MCMC_EXCHANGEUPDATER_H


#include "kpl_mcmc_dirichlet.h"


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class ExchangeabilityUpdater : public DirichletUpdater {

public:
  typedef std::shared_ptr< ExchangeabilityUpdater > SharedPtr;

  ExchangeabilityUpdater(QMatrix::SharedPtr qmatrix);
  ~ExchangeabilityUpdater();

private:

  void                            pullFromModel();
  void                            pushToModel();

  QMatrix::SharedPtr              _qmatrix;

};


} // phylogenetic
} // kellerberrin

#endif //KPL_MCMC_EXCHANGEUPDATER_H
