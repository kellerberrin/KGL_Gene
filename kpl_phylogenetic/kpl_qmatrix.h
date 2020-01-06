//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_QMATRIX_H
#define KPL_QMATRIX_H

#include "kpl_geneticcode.h"
#include "kpl_xstrom.h"

#include <Eigen/Dense>

#include <algorithm>
#include <vector>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class QMatrix {

public:
  typedef std::vector<double>             freq_xchg_t;
  typedef std::shared_ptr<freq_xchg_t>    freq_xchg_ptr_t;
  typedef double                          omega_t;
  typedef std::shared_ptr<omega_t>        omega_ptr_t;
  typedef boost::shared_ptr<QMatrix>      SharedPtr;

  QMatrix();
  virtual                                 ~QMatrix();

  virtual void                            clear() = 0;

  virtual void                            setEqualStateFreqs(freq_xchg_ptr_t freq_ptr) = 0;
  virtual void                            setStateFreqsSharedPtr(freq_xchg_ptr_t freq_ptr) = 0;
  virtual void                            setStateFreqs(freq_xchg_t & freq) = 0;
  virtual freq_xchg_ptr_t                 getStateFreqsSharedPtr() = 0;
  virtual const double *                  getStateFreqs() const = 0;
  void                                    fixStateFreqs(bool is_fixed);
  bool                                    isFixedStateFreqs() const;

  virtual void                            setEqualExchangeabilities(freq_xchg_ptr_t xchg_ptr) = 0;
  virtual void                            setExchangeabilitiesSharedPtr(freq_xchg_ptr_t xchg) = 0;
  virtual void                            setExchangeabilities(freq_xchg_t & xchg) = 0;
  virtual freq_xchg_ptr_t                 getExchangeabilitiesSharedPtr() = 0;
  virtual const double *                  getExchangeabilities() const = 0;
  void                                    fixExchangeabilities(bool is_fixed);
  bool                                    isFixedExchangeabilities() const;

  virtual void                            setOmegaSharedPtr(omega_ptr_t omega) = 0;
  virtual void                            setOmega(omega_t omega) = 0;
  virtual omega_ptr_t                     getOmegaSharedPtr() = 0;
  virtual double                          getOmega() const = 0;
  void                                    fixOmega(bool is_fixed);
  bool                                    isFixedOmega() const;

  virtual const double *                  getEigenvectors() const = 0;
  virtual const double *                  getInverseEigenvectors() const = 0;
  virtual const double *                  getEigenvalues() const = 0;

  void                                    setActive(bool activate);

protected:

  virtual void                            recalcRateMatrix() = 0;
  void                                    normalizeFreqsOrExchangeabilities(freq_xchg_ptr_t v);

  bool                                    _is_active;
  bool                                    _state_freqs_fixed;
  bool                                    _exchangeabilities_fixed;
  bool                                    _omega_fixed;
};


} // phylogenetic
} // kellerberrin


#endif // KPL_QMATRIX_H
