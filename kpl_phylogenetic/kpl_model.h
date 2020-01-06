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


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Likelihood;

class Model {

//  friend class Likelihood;

public:

  using asrv_vect_t = std::vector<ASRV::SharedPtr>;
  using qmat_vect_t = std::vector<QMatrix::SharedPtr>;
  using subset_sizes_t = std::vector<unsigned>;
  using subset_datatype_t = std::vector<DataType>;
  using subset_relrate_vect_t = std::vector<double>;
  using SharedPtr = boost::shared_ptr<Model>;
  using state_freq_params_t = std::vector<QMatrix::SharedPtr>;
  using exchangeability_params_t = std::vector<QMatrix::SharedPtr>;
  using omega_params_t = std::vector<QMatrix::SharedPtr>;
  using ratevar_params_t = std::vector<ASRV::SharedPtr>;
  using pinvar_params_t = std::vector<ASRV::SharedPtr>;

  Model();
  ~Model();

  void activate();
  void inactivate();

  [[nodiscard]] std::string describeModel();

  void setSubsetDataTypes(const subset_datatype_t & datatype_vect);

  void setSubsetRateVar(ASRV::ratevar_ptr_t ratevar, unsigned subset, bool fixed);
  void setSubsetPinvar(ASRV::pinvar_ptr_t pinvar, unsigned subset, bool fixed);
  void setSubsetExchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed);
  void setSubsetStateFreqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed);
  void setSubsetOmega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed);

  void setSubsetRelRates(subset_relrate_vect_t & relrates, bool fixed);
  subset_relrate_vect_t& getSubsetRelRates();
  [[nodiscard]] bool isFixedSubsetRelRates() const;
  [[nodiscard]] double calcNormalizingConstantForSubsetRelRates() const;

  void setTreeIndex(unsigned i, bool fixed);
  [[nodiscard]] unsigned getTreeIndex() const;
  [[nodiscard]] bool isFixedTree() const;

  void setTopologyPriorOptions(bool allow_polytomies, bool resclass, double C);
  [[nodiscard]] bool isResolutionClassTopologyPrior() const;
  [[nodiscard]] double getTopologyPriorC() const;
  [[nodiscard]] bool isAllowPolytomies() const;

  [[nodiscard]] unsigned getNumSubsets() const;
  [[nodiscard]] unsigned getNumSites() const;
  [[nodiscard]] unsigned getSubsetNumSites(unsigned subset) const;
  [[nodiscard]] const QMatrix& getQMatrix(unsigned subset) const;
  [[nodiscard]] const ASRV& getASRV(unsigned subset) const;

  void setSubsetIsInvarModel(bool is_invar, unsigned subset);
  [[nodiscard]] bool getSubsetIsInvarModel(unsigned subset) const;

  void setSubsetSizes(const subset_sizes_t nsites_vect);
  subset_sizes_t& getSubsetSizes();

  void setSubsetNumPatterns(const subset_sizes_t npatterns_vect);
  [[nodiscard]] unsigned getSubsetNumPatterns(unsigned subset) const;

  void setSubsetNumCateg(unsigned ncateg, unsigned subset);
  [[nodiscard]] unsigned getSubsetNumCateg(unsigned subset) const;

  [[nodiscard]] state_freq_params_t& getStateFreqParams();
  [[nodiscard]] exchangeability_params_t& getExchangeabilityParams();
  [[nodiscard]] omega_params_t& getOmegaParams();
  [[nodiscard]] ratevar_params_t& getRateVarParams();
  [[nodiscard]] pinvar_params_t& getPinvarParams();

  int setBeagleEigenDecomposition(int beagle_instance, unsigned subset, unsigned instance_subset);
  int setBeagleStateFrequencies(int beagle_instance, unsigned subset, unsigned instance_subset);
  int setBeagleAmongSiteRateVariationRates(int beagle_instance, unsigned subset, unsigned instance_subset);
  int setBeagleAmongSiteRateVariationProbs(int beagle_instance, unsigned subset, unsigned instance_subset);


  [[nodiscard]] std::string paramNamesAsString(std::string sep) const;
  [[nodiscard]] std::string paramValuesAsString(std::string sep) const;


private:

  void                        clear();

  unsigned                    _num_subsets;
  unsigned                    _num_sites;
  subset_sizes_t              _subset_sizes;
  subset_sizes_t              _subset_npatterns;
  subset_datatype_t           _subset_datatypes;
  qmat_vect_t                 _qmatrix;
  asrv_vect_t                 _asrv;

  bool                        _tree_index;
  bool                        _tree_fixed;

  bool                        _allow_polytomies;
  bool                        _resolution_class_prior;
  double                      _topo_prior_C;

  bool                        _subset_relrates_fixed;
  subset_relrate_vect_t       _subset_relrates;

  state_freq_params_t         _state_freq_params;
  exchangeability_params_t    _exchangeability_params;
  omega_params_t              _omega_params;
  ratevar_params_t            _ratevar_params;
  pinvar_params_t             _pinvar_params;

};


} // phylogenetic
} // kellerberrin


#endif //KPL_MODEL_H
