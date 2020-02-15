//
// Created by kellerberrin on 20/12/19.
//

#ifndef KPL_MCMC_POLYTOMYPRIOR_H
#define KPL_MCMC_POLYTOMYPRIOR_H


#include "kpl_tree.h"
#include "kpl_random.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


class PolytomyTopoPriorCalculator {

public:

  typedef std::shared_ptr<PolytomyTopoPriorCalculator> SharedPtr;

  PolytomyTopoPriorCalculator();
  ~PolytomyTopoPriorCalculator();

  bool isResolutionClassPrior() const;
  bool isPolytomyPrior() const;

  bool isRooted() const;
  bool isUnrooted() const;

  void setNTax(unsigned n);
  unsigned getNTax() const;

  void chooseRooted();
  void chooseUnrooted();

  double getLogCount(unsigned n, unsigned m);
  double getLogSaturatedCount(unsigned n);
  double getLogTotalCount(unsigned n);
  std::vector<double> getLogCounts();

  std::vector<double> getCountsVect();
  std::vector<int> getNFactorsVect();

  void chooseResolutionClassPrior();
  void choosePolytomyPrior();

  void setC(double c);
  double getC() const;

  void setLogScalingFactor(double lnf);
  double getLogScalingFactor() const;

  virtual double getLogTopoProb(Tree::SharedPtr t);

  double getLogTopologyPrior(unsigned m);
  double getLogNormalizedTopologyPrior(unsigned m);
  double getLogNormConstant();

  std::vector<double> getTopoPriorVect();
  std::vector<double> getRealizedResClassPriorsVect();

  unsigned sample(Lot::SharedPtr rng);

  void reset();

private:

  void recalcCountsAndPriorsImpl(unsigned n);
  void recalcPriorsImpl();

  unsigned _ntax;
  bool _is_rooted;
  bool _is_resolution_class_prior;
  double _C;

  bool _topo_priors_dirty;
  bool _counts_dirty;

  double _log_scaling_factor;
  std::vector<int> _nfactors;
  std::vector<double> _counts;
  double _log_total_count;
  std::vector<double> _topology_prior;
};


} // end namespace


#endif // KPL_MCMC_POLYTOMYPRIOR_H
