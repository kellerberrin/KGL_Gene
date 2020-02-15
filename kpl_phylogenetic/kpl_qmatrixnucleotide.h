//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_QMATRIXNUCLEOTIDE_H
#define KPL_QMATRIXNUCLEOTIDE_H


#include "kpl_qmatrix.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


class QMatrixNucleotide : public QMatrix {

public:
  typedef Eigen::Matrix<double, 4, 4, Eigen::RowMajor>    eigenMatrix4d_t;
  typedef Eigen::Vector4d                                 eigenVector4d_t;

  QMatrixNucleotide();
  ~QMatrixNucleotide();

  void                        clear();

  void                        setEqualStateFreqs(freq_xchg_ptr_t freq_ptr);
  void                        setStateFreqsSharedPtr(freq_xchg_ptr_t freq_ptr);
  void                        setStateFreqs(freq_xchg_t & freqs);
  freq_xchg_ptr_t             getStateFreqsSharedPtr();
  const double *              getStateFreqs() const;

  void                        setEqualExchangeabilities(freq_xchg_ptr_t xchg_ptr);
  void                        setExchangeabilitiesSharedPtr(freq_xchg_ptr_t xchg_ptr);
  void                        setExchangeabilities(freq_xchg_t & xchg);
  freq_xchg_ptr_t             getExchangeabilitiesSharedPtr();
  const double *              getExchangeabilities() const;

  void                        setOmegaSharedPtr(omega_ptr_t omega_ptr);
  void                        setOmega(omega_t omega);
  omega_ptr_t                 getOmegaSharedPtr();
  double                      getOmega() const;

  const double *              getEigenvectors() const;
  const double *              getInverseEigenvectors() const;
  const double *              getEigenvalues() const;


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:

  virtual void                recalcRateMatrix();

private:

  // workspaces for computing eigenvectors/eigenvalues
  eigenMatrix4d_t             _sqrtPi;
  eigenMatrix4d_t             _sqrtPiInv;
  eigenMatrix4d_t             _Q;
  eigenMatrix4d_t             _eigenvectors;
  eigenMatrix4d_t             _inverse_eigenvectors;
  eigenVector4d_t             _eigenvalues;

  freq_xchg_ptr_t             _state_freqs;
  freq_xchg_ptr_t             _exchangeabilities;
};


} // end namespace



#endif // KPL_QMATRIXNUCLEOTIDE_H
