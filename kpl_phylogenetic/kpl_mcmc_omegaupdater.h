//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_OMEGAUPDATER_H
#define KPL_MCMC_OMEGAUPDATER_H

#include "kpl_model.h"
#include "kpl_mcmc_updater.h"
#include "kpl_qmatrix.h"



namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class OmegaUpdater : public Updater {

public:
  typedef std::shared_ptr<OmegaUpdater> SharedPtr;

  explicit OmegaUpdater(QMatrix::SharedPtr q);
  ~OmegaUpdater() override;

  void clear() override;
  double getCurrentPoint() const;

  // mandatory overrides of pure virtual functions
  double calcLogPrior() override;
  void revert() override;
  void proposeNewState() override;

private:

  double _prev_point;
  QMatrix::SharedPtr _q;
};


} // phylogenetic
} // kellerberrin




#endif //KPL_MCMC_OMEGAUPDATER_H
