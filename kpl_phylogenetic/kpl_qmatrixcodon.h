//
// Created by kellerberrin on 13/12/19.
//

#ifndef KPL_QMATRIXCODON_H
#define KPL_QMATRIXCODON_H



#include "kpl_qmatrix.h"


namespace kellerberrin::phylogenetic {   //  organization level namespace


class QMatrixCodon : public QMatrix {

public:
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>      eigenMatrixXd_t;
  typedef Eigen::VectorXd                                                             eigenVectorXd_t;

  QMatrixCodon(GeneticCode::SharedPtr gcode);
  ~QMatrixCodon();

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

protected:

  virtual void                recalcRateMatrix();

private:

  // workspaces for computing eigenvectors/eigenvalues
  eigenMatrixXd_t             _sqrtPi;
  eigenMatrixXd_t             _sqrtPiInv;
  eigenMatrixXd_t             _Q;
  eigenMatrixXd_t             _eigenvectors;
  eigenMatrixXd_t             _inverse_eigenvectors;
  eigenVectorXd_t             _eigenvalues;

  freq_xchg_ptr_t             _state_freqs;
  omega_ptr_t                 _omega;

  std::vector<std::string>    _codons;
  std::vector<unsigned>       _amino_acids;

  GeneticCode::SharedPtr      _genetic_code;
};


} // end namespace


#endif //KPL_QMATRIXCODON_H
