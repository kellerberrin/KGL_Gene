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

  explicit GammaRateVarUpdater(ASRV::SharedPtr asrv);
  ~GammaRateVarUpdater() override;

  void clear() override;
  double  getCurrentPoint() const;

  // mandatory overrides of pure virtual functions
  double calcLogPrior() override;
  void revert() override;
  void proposeNewState() override;

private:

  double _prev_point;
  ASRV::SharedPtr _asrv;

};


} // phylogenetic
} // kellerberrin



#endif // KPL_MCMC_GAMMA_H
