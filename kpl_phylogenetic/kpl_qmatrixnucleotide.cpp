//
// Created by kellerberrin on 13/12/19.
//

#include "kpl_qmatrixnucleotide.h"

#include  <numeric>

namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here


kpl::QMatrixNucleotide::QMatrixNucleotide() {

  //std::cout << "Constructing a QMatrixNucleotide object" << std::endl;
  clear();

}


kpl::QMatrixNucleotide::~QMatrixNucleotide() {

  //std::cout << "Destroying a QMatrixNucleotide object" << std::endl;

}


void kpl::QMatrixNucleotide::clear() {

  QMatrix::clear();

  QMatrix::freq_xchg_t xchg = {1,1,1,1,1,1};
  _exchangeabilities = std::make_shared<QMatrix::freq_xchg_t>(xchg);

  QMatrix::freq_xchg_t freq_vect = {0.25, 0.25, 0.25, 0.25};
  _state_freqs = std::make_shared<QMatrix::freq_xchg_t>(freq_vect);

  recalcRateMatrix();

}


kpl::QMatrix::freq_xchg_ptr_t kpl::QMatrixNucleotide::getExchangeabilitiesSharedPtr() {

  return _exchangeabilities;

}


kpl::QMatrix::freq_xchg_ptr_t kpl::QMatrixNucleotide::getStateFreqsSharedPtr() {

  return _state_freqs;

}


kpl::QMatrix::omega_ptr_t kpl::QMatrixNucleotide::getOmegaSharedPtr() {

  assert(false);
  return nullptr;

}


const double* kpl::QMatrixNucleotide::getEigenvectors() const {

  return _eigenvectors.data();

}


const double* kpl::QMatrixNucleotide::getInverseEigenvectors() const {

  return _inverse_eigenvectors.data();

}


const double* kpl::QMatrixNucleotide::getEigenvalues() const {

  return _eigenvalues.data();

}


const double* kpl::QMatrixNucleotide::getExchangeabilities() const {

  return &(*_exchangeabilities)[0];

}


const double* kpl::QMatrixNucleotide::getStateFreqs() const {

  return &(*_state_freqs)[0];

}


double kpl::QMatrixNucleotide::getOmega() const {

  assert(false);
  return 0.0;

}


void kpl::QMatrixNucleotide::setEqualExchangeabilities(QMatrix::freq_xchg_ptr_t xchg_ptr) {

  _exchangeabilities = xchg_ptr;
  _exchangeabilities->assign(6, 1.0/6.0);
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setExchangeabilitiesSharedPtr(QMatrix::freq_xchg_ptr_t xchg_ptr) {

  if (xchg_ptr->size() != 6) {

    throw XStrom(boost::format("Expecting 6 exchangeabilities and got %d: perhaps you meant to specify a subset data type other than nucleotide") % xchg_ptr->size());

  }

  _exchangeabilities = xchg_ptr;
  normalizeFreqsOrExchangeabilities(_exchangeabilities);
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setExchangeabilities(QMatrix::freq_xchg_t & xchg) {

  if (xchg.size() != 6) {

    throw XStrom(boost::format("Expecting 6 exchangeabilities and got %d: perhaps you meant to specify a subset data type other than nucleotide") % xchg.size());

  }

  std::copy(xchg.begin(), xchg.end(), _exchangeabilities->begin());
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setEqualStateFreqs(QMatrix::freq_xchg_ptr_t freq_ptr) {

  _state_freqs = freq_ptr;
  _state_freqs->assign(4, 0.25);
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setStateFreqsSharedPtr(QMatrix::freq_xchg_ptr_t freq_ptr) {

  if (freq_ptr->size() != 4) {

    throw XStrom(boost::format("Expecting 4 state frequencies and got %d: perhaps you meant to specify a subset data type other than nucleotide") % freq_ptr->size());

  }

  double sum_of_freqs = std::accumulate(freq_ptr->begin(), freq_ptr->end(), 0.0);

  if (std::fabs(sum_of_freqs - 1.0) > 0.001) {

    throw XStrom(boost::format("Expecting sum of 4 state frequencies to be 1, but instead got %g") % sum_of_freqs);

  }

  _state_freqs = freq_ptr;
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setStateFreqs(QMatrix::freq_xchg_t & freqs) {

  if (freqs.size() != 4) {

    throw XStrom(boost::format("Expecting 4 state frequencies and got %d: perhaps you meant to specify a subset data type other than nucleotide") % freqs.size());

  }

  std::copy(freqs.begin(), freqs.end(), _state_freqs->begin());
  recalcRateMatrix();

}


void kpl::QMatrixNucleotide::setOmegaSharedPtr(QMatrix::omega_ptr_t omega_ptr) {

  assert(false);

}


void kpl::QMatrixNucleotide::setOmega(QMatrix::omega_t omega) {

  assert(false);

}


void kpl::QMatrixNucleotide::recalcRateMatrix() {

  // Must have assigned both _state_freqs and _exchangeabilities to recalculate rate matrix
  if (!_is_active || !(_state_freqs && _exchangeabilities)) {

    return;

  }

  double piA = (*_state_freqs)[0];
  double piC = (*_state_freqs)[1];
  double piG = (*_state_freqs)[2];
  double piT = (*_state_freqs)[3];

  Eigen::Map<const Eigen::Array4d> tmp(_state_freqs->data());
  _sqrtPi = tmp.sqrt().matrix().asDiagonal();
  _sqrtPiInv = _sqrtPi.inverse();

  assert(_exchangeabilities->size() == 6);

  double rAC = (*_exchangeabilities)[0];
  double rAG = (*_exchangeabilities)[1];
  double rAT = (*_exchangeabilities)[2];
  double rCG = (*_exchangeabilities)[3];
  double rCT = (*_exchangeabilities)[4];
  double rGT = (*_exchangeabilities)[5];

  double inverse_scaling_factor = piA*(rAC*piC + rAG*piG + rAT*piT) + piC*(rAC*piA + rCG*piG + rCT*piT) + piG*(rAG*piA + rCG*piC + rGT*piT) + piT*(rAT*piA + rCT*piC + rGT*piG);
  double scaling_factor = 1.0/inverse_scaling_factor;

  _Q(0,0) = -scaling_factor*(rAC*piC + rAG*piG + rAT*piT);
  _Q(0,1) = scaling_factor*rAC*piC;
  _Q(0,2) = scaling_factor*rAG*piG;
  _Q(0,3) = scaling_factor*rAT*piT;

  _Q(1,0) = scaling_factor*rAC*piA;
  _Q(1,1) = -scaling_factor*(rAC*piA + rCG*piG + rCT*piT);
  _Q(1,2) = scaling_factor*rCG*piG;
  _Q(1,3) = scaling_factor*rCT*piT;

  _Q(2,0) = scaling_factor*rAG*piA;
  _Q(2,1) = scaling_factor*rCG*piC;
  _Q(2,2) = -scaling_factor*(rAG*piA + rCG*piC + rGT*piT);
  _Q(2,3) = scaling_factor*rGT*piT;

  _Q(3,0) = scaling_factor*rAT*piA;
  _Q(3,1) = scaling_factor*rCT*piC;
  _Q(3,2) = scaling_factor*rGT*piG;
  _Q(3,3) = -scaling_factor*(rAT*piA + rCT*piC + rGT*piG);

  // S is a symmetric matrix
  eigenMatrix4d_t S = eigenMatrix4d_t(_sqrtPi*_Q*_sqrtPiInv);

  // Can use efficient eigensystem solver because S is symmetric
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(S);

  if (solver.info() != Eigen::Success) {

    throw XStrom("Error in the calculation of eigenvectors and eigenvalues of the GTR rate matrix");

  }

  _eigenvectors           = _sqrtPiInv*solver.eigenvectors();
  _inverse_eigenvectors   = solver.eigenvectors().transpose()*_sqrtPi;
  _eigenvalues            = solver.eigenvalues();

}