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

  OmegaUpdater(QMatrix::SharedPtr q);
  ~OmegaUpdater();

  virtual void                clear();
  double                      getCurrentPoint() const;

  // mandatory overrides of pure virtual functions
  virtual double              calcLogPrior();
  virtual void                revert();
  virtual void                proposeNewState();

private:

  double                      _prev_point;
  QMatrix::SharedPtr          _q;
};


} // phylogenetic
} // kellerberrin




#endif //KPL_MCMC_OMEGAUPDATER_H
