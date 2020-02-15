//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_MCMC_GAMMA_H
#define KPL_MCMC_GAMMA_H


#include "kpl_model.h"
#include "kpl_mcmc_updater.h"
#include "kpl_asrv.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


class GammaRateVarUpdater : public Updater {

public:

  typedef std::shared_ptr<GammaRateVarUpdater> SharedPtr;

  explicit GammaRateVarUpdater(std::shared_ptr<ASRV> asrv);
  ~GammaRateVarUpdater() override;

  void clear() override;
  double  getCurrentPoint() const;

  // mandatory overrides of pure virtual functions
  double calcLogPrior() override;
  void revert() override;
  void proposeNewState() override;

private:

  double _prev_point;
  std::shared_ptr<ASRV> _asrv;

};


} // end namespace


#endif // KPL_MCMC_GAMMA_H
