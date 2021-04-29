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


namespace kellerberrin::phylogenetic {   //  organization level namespace

class Chain;

class Updater {


public:
  typedef std::shared_ptr<Updater>    SharedPtr;


  Updater();
  virtual ~Updater() = default;

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
  void name(const std::string& update) { name_ = update; }
  void logHastingsRatio(double log_hastings_ratio) { log_hastings_ratio_ = log_hastings_ratio; }
  void logJacobian(double log_jacobian) { log_jacobian_ = log_jacobian; }
  void lambda(double lambda) { lambda_ = lambda; }


  // Access
  [[nodiscard]] double getLambda() const;
  [[nodiscard]] double getWeight() const;
  [[nodiscard]] double getAcceptPct() const;
  [[nodiscard]] double getNumUpdates() const;
  [[nodiscard]] std::string getUpdaterName() const;
  [[nodiscard]] const std::string& name() const { return name_; }
  [[nodiscard]] double probability() const { return prob_; }
  [[nodiscard]] const std::vector<double>& priorParameters() { return prior_parameters_; }
  [[nodiscard]] static double logZero() { return LOG_ZERO_; }
  [[nodiscard]] double lambda() const { return lambda_; }
  [[nodiscard]] double logHastingsRatio() const { return log_hastings_ratio_; }
  [[nodiscard]] double logJacobian() const { return log_jacobian_; }


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
  [[nodiscard]] const std::shared_ptr<Lot>& lot() const { return lot_ptr_; }
  [[nodiscard]] const std::shared_ptr<TreeManip> treeManipulator() const { return tree_manipulator_ptr_; }

private:

  virtual void tune(bool accepted);
  virtual void revert() = 0;
  virtual void proposeNewState() = 0;

  std::shared_ptr<Lot> lot_ptr_;
  std::shared_ptr<Likelihood> likelihood_ptr_;
  std::shared_ptr<TreeManip> tree_manipulator_ptr_;
  std::string name_;
  double weight_;
  double prob_;
  double lambda_;
  double log_hastings_ratio_;
  double log_jacobian_;
  double target_acceptance_;
  unsigned n_accepts_;
  unsigned n_attempts_;
  bool tuning_;
  std::vector<double> prior_parameters_;

  double heating_power_;

  mutable PolytomyTopoPriorCalculator topo_prior_calculator_;

// The largest negative value.
  static constexpr double const LOG_ZERO_ = std::numeric_limits<double>::lowest();

};


} // end namespace


#endif // KPL_MCMC_UPDATER_H
