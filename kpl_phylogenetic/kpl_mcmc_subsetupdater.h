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

  explicit SubsetRelRateUpdater(const std::shared_ptr<Model>& model_ptr);
  ~SubsetRelRateUpdater() override = default;

  double calcLogPrior() override;

private:

  void pullFromModel() override;
  void pushToModel() override;

  std::shared_ptr<Model> model_ptr_;

};



} // end namespace


#endif // KPL_MCMC_SUBSETUPDATER_H
