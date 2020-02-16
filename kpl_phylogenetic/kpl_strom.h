//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_STROM_H
#define KPL_STROM_H

#include "kel_exec_env.h"

#include "kpl_geneticdata.h"
#include "kpl_treesummary.h"
#include "kpl_partition.h"
#include "kpl_likelihood.h"
#include "kpl_model.h"
#include "kpl_random.h"
#include "kpl_mcmc_chain.h"
#include "kpl_mcmc_output.h"


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>


namespace kellerberrin::phylogenetic {   //  organization level namespace



class Strom {

public:
  Strom() { clear(); }
  ~Strom() = default;

  void parseCommandLine(int argc, const char **argv);
  void executeApp();

  static constexpr const char* VERSION = "0.1";
  static constexpr const char* MODULE_NAME = "kpl_phyloTree";

private:

  void clear();
  bool processAssignmentString(std::shared_ptr<Model> model_ptr, const std::string & which, const std::string & definition);
  void handleAssignmentStrings(std::shared_ptr<Model> model_ptr, const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition);
  bool splitAssignmentString(const std::string & definition, std::vector<std::string> & vector_of_subset_names, std::vector<double>  & vector_of_values);
  void sample(unsigned iter, Chain & chain);

  void readData();
  void readTrees();
  void showPartitionInfo();
  void showBeagleInfo();
  void showMCMCInfo();
  void calcHeatingPowers();
  void initChains();
  void startTuningChains();
  void stopTuningChains();
  void stepChains(unsigned iteration, bool sampling);
  void swapChains();
  void stopChains();
  void swapSummary() const;
  void showChainTuningInfo() const;

  double _expected_log_likelihood;

  double _topo_prior_C;
  bool _allow_polytomies;
  bool _resolution_class_prior;

  std::string _data_file_name;
  std::string _tree_file_name;
  Partition::SharedPtr _partition;

  Data::SharedPtr _data;
  std::vector<Likelihood::SharedPtr> _likelihoods;
  TreeSummary::SharedPtr _tree_summary;
  Lot::SharedPtr _lot;

  unsigned _random_seed;
  unsigned _num_iter;
  unsigned _print_freq;
  unsigned _sample_freq;

  unsigned _num_burnin_iter;
  unsigned _num_chains;
  double _heating_lambda;
  bool _using_stored_data;
  std::vector<Chain> _chains;
  std::vector<double> _heating_powers;
  std::vector<unsigned> _swaps;

  bool _use_gpu;
  bool _ambig_missing;
  bool _use_underflow_scaling;

  OutputManager::SharedPtr _output_manager;


};



} // end namespace


#endif //KPL_STROM_H
