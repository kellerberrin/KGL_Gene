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


namespace kellerberrin::phylogenetic {   //  organization::project level namespace


class Chain {


public:

  using UpdaterVector = std::vector<std::shared_ptr<Updater>>;
  using SharedPtr = std::shared_ptr<Chain>;

  Chain();
  ~Chain() = default;

  void startTuning();
  void stopTuning();

  void setTreeFromNewick(const std::string& newick);
  [[nodiscard]] size_t createUpdaters(const std::shared_ptr<Model>& model_ptr,
                                      const std::shared_ptr<Lot>& lot,
                                      const std::shared_ptr<Likelihood>& likelihood);

  [[nodiscard]] const std::shared_ptr<TreeManip>& getTreeManip() const { return tree_manip_ptr_; }
  [[nodiscard]] const std::shared_ptr<Model>& getModel() const { return model_ptr_; }
  [[nodiscard]] double getLogLikelihood() const { return log_likelihood_; }


  void setHeatingPower(double p);
  [[nodiscard]] double getHeatingPower() const { return heating_power_; }

  void setChainIndex(unsigned idx);
  [[nodiscard]] double getChainIndex() const;

  [[nodiscard]] std::vector<std::string> getUpdaterNames() const;
  [[nodiscard]] std::vector<double> getAcceptPercentages() const;
  [[nodiscard]] std::vector<unsigned> getNumUpdates() const;
  [[nodiscard]] std::vector<double> getLambdas() const;
  void setLambdas(std::vector<double> & v);

  [[nodiscard]] double calcLogLikelihood() const;
  [[nodiscard]] double calcLogJointPrior() const;

  void start();
  void stop();
  void nextStep(int iteration);

private:

  std::shared_ptr<Model> model_ptr_;
  std::shared_ptr<Lot> lot_ptr_;
  std::shared_ptr<TreeManip> tree_manip_ptr_;

  UpdaterVector updaters_;

  size_t chain_index_{0};
  double heating_power_{1.0};
  double log_likelihood_{0.0};

};


} // end namespace


#endif //KPL_MCMC_CHAIN_H
