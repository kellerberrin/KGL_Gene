//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_STATEFREQUPDATER_H
#define KPL_MCMC_STATEFREQUPDATER_H


#include "kpl_mcmc_dirichlet.h"


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class StateFreqUpdater : public DirichletUpdater {

public:

  typedef std::shared_ptr< StateFreqUpdater > SharedPtr;

  StateFreqUpdater(QMatrix::SharedPtr qmatrix);
  ~StateFreqUpdater();

private:

  void                            pullFromModel();
  void                            pushToModel();

  QMatrix::SharedPtr              _qmatrix;

};


} // phylogenetic
} // kellerberrin




#endif //KPL_MCMC_STATEFREQUPDATER_H
