//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_MCMC_CHAIN_H
#define KPL_MCMC_CHAIN_H

#include "kpl_random.h"
#include "kpl_geneticdata.h"
#include "kpl_tree.h"
#include "kpl_model.h"
#include "kpl_likelihood.h"
#include "kpl_treemanip.h"
#include "kpl_mcmc_updater.h"
#include "kpl_mcmc_gamma.h"
#include "kpl_mcmc_treeupdater.h"
#include "kpl_mcmc_treelength.h"
#include "kpl_mcmc_polytomyupdater.h"



#include <boost/format.hpp>

#include <memory>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Chain {

//  friend class Likelihood;

public:

  typedef std::vector<Updater::SharedPtr> updater_vect_t;
  typedef std::shared_ptr<Chain>          SharedPtr;

  Chain();
  ~Chain();

  void                                    clear();

  void                                    startTuning();
  void                                    stopTuning();

  void                                    setTreeFromNewick(std::string & newick);
  unsigned                                createUpdaters(Model::SharedPtr model, Lot::SharedPtr lot, Likelihood::SharedPtr likelihood);

  TreeManip::SharedPtr                    getTreeManip();
  Model::SharedPtr                        getModel();
  double                                  getLogLikelihood() const;


  void                                    setHeatingPower(double p);
  double                                  getHeatingPower() const;

  void                                    setChainIndex(unsigned idx);
  double                                  getChainIndex() const;

  Updater::SharedPtr                      findUpdaterByName(std::string name);
  std::vector<std::string>                getUpdaterNames() const;
  std::vector<double>                     getAcceptPercentages() const;
  std::vector<unsigned>                   getNumUpdates() const;
  std::vector<double>                     getLambdas() const;
  void                                    setLambdas(std::vector<double> & v);

  double                                  calcLogLikelihood() const;
  double                                  calcLogJointPrior() const;

  void                                    start();
  void                                    stop();
  void                                    nextStep(int iteration);

private:

  Model::SharedPtr                        _model;
  Lot::SharedPtr                          _lot;
  TreeManip::SharedPtr                    _tree_manipulator;

  updater_vect_t                          _updaters;

  unsigned                                _chain_index;
  double                                  _heating_power;
  double                                  _log_likelihood;
};


} // phylogenetic
} // kellerberrin


#endif //KPL_MCMC_CHAIN_H
