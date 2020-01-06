//
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_SUBSETUPDATER_H
#define KPL_MCMC_SUBSETUPDATER_H



#include "kpl_mcmc_dirichlet.h"



namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class SubsetRelRateUpdater : public DirichletUpdater {
public:
  typedef std::shared_ptr< SubsetRelRateUpdater > SharedPtr;

  SubsetRelRateUpdater(Model::SharedPtr model);
  ~SubsetRelRateUpdater();

  double                          calcLogPrior();

private:

  void                            pullFromModel();
  void                            pushToModel();

  Model::SharedPtr                _model;

};



} // phylogenetic
} // kellerberrin


#endif // KPL_MCMC_SUBSETUPDATER_H
