//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_MODEL_H
#define KPL_MODEL_H


#include "kpl_genetictype.h"
#include "kpl_qmatrix.h"
#include "kpl_asrv.h"

#include "libhmsbeagle/beagle.h"
#include <boost/format.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <Eigen/Dense>

#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>


namespace kellerberrin::phylogenetic {   //  organization::project level namespace


class Likelihood;

class Model {

//  friend class Likelihood;

public:

  using asrv_vect_t = std::vector<std::shared_ptr<ASRV> >;
  using ratevar_params_t = std::vector<std::shared_ptr<ASRV> >;
  using pinvar_params_t = std::vector<std::shared_ptr<ASRV> >;

  using qmat_vect_t = std::vector<QMatrix::SharedPtr>;
  using subset_sizes_t = std::vector<unsigned>;
  using subset_datatype_t = std::vector<DataType>;
  using subset_relrate_vect_t = std::vector<double>;
  using state_freq_params_t = std::vector<QMatrix::SharedPtr>;
  using exchangeability_params_t = std::vector<QMatrix::SharedPtr>;
  using omega_params_t = std::vector<QMatrix::SharedPtr>;

  Model() { clear(); }
  ~Model() = default;

  void activate();
  void inactivate();

// Modify
  void setSubsetDataTypes(const subset_datatype_t & datatype_vect);
  void setSubsetRateVar(std::shared_ptr<double> ratevar_ptr, unsigned subset, bool fixed);
  void setSubsetPinvar(std::shared_ptr<double> pinvar_ptr, unsigned subset, bool fixed);
  void setSubsetExchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed);
  void setSubsetStateFreqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed);
  void setSubsetOmega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed);
  void setTreeIndex(unsigned i, bool fixed);
  void setSubsetRelRates(subset_relrate_vect_t & relrates, bool fixed);
  void setTopologyPriorOptions(bool allow_polytomies, bool resclass, double C);
  void setSubsetIsInvarModel(bool is_invar, unsigned subset);
  void setSubsetSizes(const subset_sizes_t nsites_vect);
  void setSubsetNumPatterns(const subset_sizes_t npatterns_vect);
  void setSubsetNumCateg(unsigned ncateg, unsigned subset);

  std::string describeModel(); // Modifies the object - todo: separate model description from model modification.

// Access
  [[nodiscard]] bool isFixedSubsetRelRates() const { return _subset_relrates_fixed; }
  [[nodiscard]] unsigned getTreeIndex() const { return _tree_index; }
  [[nodiscard]] bool isFixedTree() const { return _tree_fixed; }
  [[nodiscard]] bool isResolutionClassTopologyPrior() const { return _resolution_class_prior; }
  [[nodiscard]] double getTopologyPriorC() const { return _topo_prior_C; }
  [[nodiscard]] bool isAllowPolytomies() const { return _allow_polytomies; }
  [[nodiscard]] unsigned getNumSubsets() const { return _num_subsets; }
  [[nodiscard]] unsigned getNumSites() const { return _num_sites; }
  [[nodiscard]] const subset_sizes_t& getSubsetSizes() const { return _subset_sizes; }
  [[nodiscard]] const state_freq_params_t& getStateFreqParams() const { return _state_freq_params; }
  [[nodiscard]] const exchangeability_params_t& getExchangeabilityParams() const { return _exchangeability_params; }
  [[nodiscard]] const omega_params_t& getOmegaParams() const { return _omega_params; }
  [[nodiscard]] const ratevar_params_t& getRateVarParams() const { return _ratevar_params; }
  [[nodiscard]] const pinvar_params_t& getPinvarParams() const { return _pinvar_params; }
  [[nodiscard]] unsigned getSubsetNumSites(unsigned subset) const;
  [[nodiscard]] const QMatrix& getQMatrix(unsigned subset) const;
  [[nodiscard]] const ASRV& getASRV(unsigned subset) const;
  [[nodiscard]] bool getSubsetIsInvarModel(unsigned subset) const;
  [[nodiscard]] double calcNormalizingConstantForSubsetRelRates() const;
  [[nodiscard]] unsigned getSubsetNumCateg(unsigned subset) const;
  [[nodiscard]] std::string paramNamesAsString(std::string sep) const;
  [[nodiscard]] std::string paramValuesAsString(std::string sep) const;
  [[nodiscard]] const subset_relrate_vect_t& getSubsetRelRates() const;


private:

  void clear();

  unsigned _num_subsets;
  unsigned _num_sites;
  subset_sizes_t _subset_sizes;
  subset_sizes_t _subset_npatterns;
  subset_datatype_t _subset_datatypes;
  qmat_vect_t _qmatrix;
  asrv_vect_t _asrv;

  bool _tree_index;
  bool _tree_fixed;

  bool _allow_polytomies;
  bool _resolution_class_prior;
  double _topo_prior_C;

  bool _subset_relrates_fixed;
  subset_relrate_vect_t _subset_relrates;

  state_freq_params_t _state_freq_params;
  exchangeability_params_t _exchangeability_params;
  omega_params_t _omega_params;
  ratevar_params_t _ratevar_params;
  pinvar_params_t _pinvar_params;

};


} // end namespace


#endif //KPL_MODEL_H
