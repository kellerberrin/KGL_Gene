//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_ASRV_H
#define KPL_ASRV_H

#include <boost/math/distributions/gamma.hpp>

#include <vector>
#include <memory>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class ASRV {

public:
  typedef std::vector<double>         rate_prob_t;
  typedef std::shared_ptr<double>     relrate_ptr_t;
  typedef std::shared_ptr<double>     ratevar_ptr_t;
  typedef std::shared_ptr<double>     pinvar_ptr_t;
  typedef std::shared_ptr<ASRV>       SharedPtr;

  ASRV();
  virtual                             ~ASRV();

  void                                clear();

  void                                setNumCateg(unsigned ncateg);
  unsigned                            getNumCateg() const;

  void                                setRateVarSharedPtr(ratevar_ptr_t ratevar);
  void                                setRateVar(double v);
  const ratevar_ptr_t                 getRateVarSharedPtr() const;
  double                              getRateVar() const;
  void                                fixRateVar(bool is_fixed);
  bool                                isFixedRateVar() const;

  void                                setPinvarSharedPtr(pinvar_ptr_t pinvar);
  void                                setPinvar(double p);
  const pinvar_ptr_t                  getPinvarSharedPtr() const;
  double                              getPinvar() const;
  void                                fixPinvar(bool is_fixed);
  bool                                isFixedPinvar() const;

  void                                setIsInvarModel(bool is_invar_model);
  bool                                getIsInvarModel() const;

  const double *                      getRates() const;
  const double *                      getProbs() const;

private:

  virtual void                        recalcASRV();

  unsigned                            _num_categ;
  bool                                _invar_model;

  ratevar_ptr_t                       _ratevar;
  pinvar_ptr_t                        _pinvar;

  bool                                _ratevar_fixed;
  bool                                _pinvar_fixed;

  rate_prob_t                         _rates;
  rate_prob_t                         _probs;
};


} // phylogenetic
} // kellerberrin


#endif //STROM_KPL_ASRV_H
