////
// Created by kellerberrin on 17/12/19.
//

#ifndef KPL_MCMC_DIRICHLET_H
#define KPL_MCMC_DIRICHLET_H


#include "kpl_mcmc_updater.h"


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Chain;

class DirichletUpdater : public Updater {

//  friend class Chain;

public:
  typedef std::vector<double>                 point_t;
  typedef std::shared_ptr< DirichletUpdater > SharedPtr;

  DirichletUpdater();
  virtual                             ~DirichletUpdater();

  void                                clear();
  virtual double                      calcLogPrior();

protected:

  virtual void                        pullFromModel() = 0;
  virtual void                        pushToModel() = 0;

  void                                proposeNewState();
  void                                revert();

  point_t                             _curr_point;
  point_t                             _prev_point;


};



} // phylogenetic
} // kellerberrin


#endif //KPL_MCMC_DIRICHLET_H
