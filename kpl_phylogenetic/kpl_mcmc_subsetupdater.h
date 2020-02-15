//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_SUBSETUPDATER_H
#define KPL_MCMC_SUBSETUPDATER_H



#include "kpl_mcmc_dirichlet.h"



namespace kellerberrin::phylogenetic {   //  organization::project level namespace


class SubsetRelRateUpdater : public DirichletUpdater {
public:
  typedef std::shared_ptr< SubsetRelRateUpdater > SharedPtr;

  SubsetRelRateUpdater(std::shared_ptr<Model> model);
  ~SubsetRelRateUpdater();

  double calcLogPrior() override;

private:

  void pullFromModel();
  void pushToModel();

  std::shared_ptr<Model> _model;

};



} // end namespace


#endif // KPL_MCMC_SUBSETUPDATER_H
