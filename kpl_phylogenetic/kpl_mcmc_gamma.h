//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_MCMC_GAMMA_H
#define KPL_MCMC_GAMMA_H


#include "kpl_model.h"
#include "kpl_mcmc_updater.h"
#include "kpl_asrv.h"


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class GammaRateVarUpdater : public Updater {

public:
  typedef std::shared_ptr<GammaRateVarUpdater> SharedPtr;

  GammaRateVarUpdater(ASRV::SharedPtr asrv);
  ~GammaRateVarUpdater();

  virtual void                clear();
  double                      getCurrentPoint() const;

  // mandatory overrides of pure virtual functions
  virtual double              calcLogPrior();
  virtual void                revert();
  virtual void                proposeNewState();

private:

  double                      _prev_point;
  ASRV::SharedPtr             _asrv;

};


} // phylogenetic
} // kellerberrin



#endif // KPL_MCMC_GAMMA_H
