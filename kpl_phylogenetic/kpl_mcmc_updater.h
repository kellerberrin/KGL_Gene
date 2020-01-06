//
// Created by kellerberrin on 16/12/19.
//

#ifndef KPL_MCMC_UPDATER_H
#define KPL_MCMC_UPDATER_H


#include "kpl_tree.h"
#include "kpl_treemanip.h"
#include "kpl_random.h"
#include "kpl_xstrom.h"
#include "kpl_likelihood.h"
#include "kpl_mcmc_polytomyprior.h"


#include <memory>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Chain;

class Updater {

//  friend class Chain;

public:
  typedef std::shared_ptr<Updater>    SharedPtr;

  TreeManip::SharedPtr                getTreeManip() const;

  Updater();
  virtual                             ~Updater();

  // Modify
  void setLikelihood(Likelihood::SharedPtr likelihood);
  void setTreeManip(TreeManip::SharedPtr treemanip);
  void setLot(Lot::SharedPtr lot);
  void setLambda(double lambda);
  void setHeatingPower(double p);
  void setTuning(bool on);
  void setTargetAcceptanceRate(double target);
  void setPriorParameters(const std::vector<double> & c);
  void setTopologyPriorOptions(bool resclass, double C);
  void setWeight(double w);
  void calcProb(double wsum);
  void name(const std::string& update) { _name = update; }
  void logHastingsRatio(double log_hastings_ratio) { _log_hastings_ratio = log_hastings_ratio; }
  void logJacobian(double log_jacobian) { _log_jacobian = log_jacobian; }
  void lambda(double lambda) { _lambda = lambda; }


  // Access
  [[nodiscard]] double getLambda() const;
  [[nodiscard]] double getWeight() const;
  [[nodiscard]] double getProb() const;
  [[nodiscard]] double getAcceptPct() const;
  [[nodiscard]] double getNumUpdates() const;
  [[nodiscard]] std::string getUpdaterName() const;
  [[nodiscard]] const std::string& name() const { return _name; }
  [[nodiscard]] double probability() const { return _prob; }
  [[nodiscard]] const std::vector<double>& priorParameters() { return _prior_parameters; }
  [[nodiscard]] static double logZero() { return _log_zero; }
  [[nodiscard]] double lambda() const { return _lambda; }
  [[nodiscard]] double logHastingsRatio() const { return _log_hastings_ratio; }
  [[nodiscard]] double logJacobian() const { return _log_jacobian; }


  virtual void clear();
  virtual double calcLogPrior() = 0;
  double calcLogTopologyPrior() const;
// Note that this function was named "double calcLogEdgeLengthPrior() const" in the source code that
// was with the tutorial. Some source code in kpl_mcmc_polytomyupdater.cc referred to this alternate
// function name in the derived class function "double kpl::PolytomyUpdater::calcLogPrior()" and was modified.
// to use "double calcEdgeLengthPrior() const" instead.
  double calcEdgeLengthPrior() const;
  double calcLogLikelihood() const;
  virtual double update(double prev_lnL);

  static double getLogZero();

protected:

  virtual void reset();
  [[nodiscard]] Lot::SharedPtr lot() { return _lot; }
  [[nodiscard]] TreeManip::SharedPtr treeManipulator() { return _tree_manipulator; }

private:

  virtual void tune(bool accepted);

  virtual void revert() = 0;
  virtual void proposeNewState() = 0;

  Lot::SharedPtr _lot;
  Likelihood::SharedPtr _likelihood;
  TreeManip::SharedPtr _tree_manipulator;
  std::string _name;
  double _weight;
  double _prob;
  double _lambda;
  double _log_hastings_ratio;
  double _log_jacobian;
  double _target_acceptance;
  unsigned _naccepts;
  unsigned _nattempts;
  bool _tuning;
  std::vector<double> _prior_parameters;

  double _heating_power;

  mutable PolytomyTopoPriorCalculator _topo_prior_calculator;

// The largest negative value.
  static constexpr double const _log_zero = std::numeric_limits<double>::lowest();

};


} // phylogenetic
} // kellerberrin


#endif // KPL_MCMC_UPDATER_H
