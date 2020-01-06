//
// Created by kellerberrin on 13/12/19.
//

#include "kpl_qmatrixcodon.h"

#include  <numeric>

namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here


kpl::QMatrixCodon::QMatrixCodon(GeneticCode::SharedPtr gcode) {
  //std::cout << "Constructing a QMatrixCodon object" << std::endl;
  assert(gcode);
  _genetic_code = gcode;
  clear();

}


kpl::QMatrixCodon::~QMatrixCodon() {
  //std::cout << "Destroying a QMatrixCodon object" << std::endl;
}


void kpl::QMatrixCodon::clear() {

  QMatrix::clear();

  unsigned nstates = _genetic_code->getNumNonStopCodons();
  _genetic_code->copyCodons(_codons);
  _genetic_code->copyAminoAcids(_amino_acids);

  QMatrix::omega_t omega = 0.1;
  _omega = std::make_shared<QMatrix::omega_t>(omega);

  QMatrix::freq_xchg_t freq_vect(nstates, 1./nstates);
  _state_freqs = std::make_shared<QMatrix::freq_xchg_t>(freq_vect);

  _sqrtPi.resize(nstates, nstates);
  _sqrtPiInv.resize(nstates, nstates);
  _Q.resize(nstates, nstates);
  _eigenvectors.resize(nstates, nstates);
  _inverse_eigenvectors.resize(nstates, nstates);
  _eigenvalues.resize(nstates);

  recalcRateMatrix();

}


kpl::QMatrix::freq_xchg_ptr_t kpl::QMatrixCodon::getExchangeabilitiesSharedPtr() {

  assert(false);
  return nullptr;

}


kpl::QMatrix::freq_xchg_ptr_t kpl::QMatrixCodon::getStateFreqsSharedPtr() {

  return _state_freqs;

}


kpl::QMatrix::omega_ptr_t kpl::QMatrixCodon::getOmegaSharedPtr() {

  return _omega;

}


const double* kpl::QMatrixCodon::getEigenvectors() const {

  return _eigenvectors.data();

}


const double* kpl::QMatrixCodon::getInverseEigenvectors() const {

  return _inverse_eigenvectors.data();

}


const double* kpl::QMatrixCodon::getEigenvalues() const {

  return _eigenvalues.data();

}

const double* kpl::QMatrixCodon::getExchangeabilities() const {

  assert(false);
  return 0;

}


const double* kpl::QMatrixCodon::getStateFreqs() const {

  return &(*_state_freqs)[0];

}


double kpl::QMatrixCodon::getOmega() const {

  return *_omega;

}


void kpl::QMatrixCodon::setEqualExchangeabilities(QMatrix::freq_xchg_ptr_t xchg_ptr) {

  assert(false);

}


void kpl::QMatrixCodon::setExchangeabilitiesSharedPtr(QMatrix::freq_xchg_ptr_t xchg_ptr) {

  assert(false);

}


void kpl::QMatrixCodon::setExchangeabilities(QMatrix::freq_xchg_t & xchg) {

  assert(false);

}


void kpl::QMatrixCodon::setEqualStateFreqs(QMatrix::freq_xchg_ptr_t freq_ptr) {

  _state_freqs = freq_ptr;
  unsigned nstates = _genetic_code->getNumNonStopCodons();
  _state_freqs->assign(nstates, 1./nstates);
  recalcRateMatrix();

}


void kpl::QMatrixCodon::setStateFreqsSharedPtr(QMatrix::freq_xchg_ptr_t freq_ptr) {
  unsigned nstates = _genetic_code->getNumNonStopCodons();
  if (freq_ptr->size() != nstates) {

    throw XStrom(boost::format("Expecting %d state frequencies and got %d: perhaps you meant to specify a subset data type other than codon") % nstates % freq_ptr->size());

  }

  double sum_of_freqs = std::accumulate(freq_ptr->begin(), freq_ptr->end(), 0.0);

  if (std::fabs(sum_of_freqs - 1.0) > 0.001) {

    throw XStrom(boost::format("Expecting sum of codon frequencies to be 1, but instead got %g") % sum_of_freqs);

  }

  _state_freqs = freq_ptr;

  normalizeFreqsOrExchangeabilities(_state_freqs);
  recalcRateMatrix();

}


void kpl::QMatrixCodon::setStateFreqs(QMatrix::freq_xchg_t & freqs) {

  unsigned nstates = _genetic_code->getNumNonStopCodons();

  if (freqs.size() != nstates) {

    throw XStrom(boost::format("Expecting %d state frequencies and got %d: perhaps you meant to specify a subset data type other than codon") % nstates % freqs.size());

  }
  std::copy(freqs.begin(), freqs.end(), _state_freqs->begin());
  recalcRateMatrix();

}


void kpl::QMatrixCodon::setOmegaSharedPtr(QMatrix::omega_ptr_t omega_ptr) {

  _omega = omega_ptr;
  recalcRateMatrix();

}


void kpl::QMatrixCodon::setOmega(QMatrix::omega_t omega) {

  *_omega = omega;
  recalcRateMatrix();

}


void kpl::QMatrixCodon::recalcRateMatrix() {
  // Must have assigned both _state_freqs and _omega to recalculate rate matrix
  if (!_is_active || !(_state_freqs && _omega)) {

    return;

  }

  unsigned nstates = _genetic_code->getNumNonStopCodons();
  assert(_state_freqs->size() == nstates);
  const double * pi = getStateFreqs();
  double omega = getOmega();

  Eigen::Map<const Eigen::ArrayXd> tmp(_state_freqs->data(), nstates);
  _sqrtPi = tmp.sqrt().matrix().asDiagonal();
  _sqrtPiInv = _sqrtPi.inverse();

  // Calculate (unscaled) instantaneous rate matrix
  _Q = Eigen::MatrixXd::Zero(nstates,nstates);

  for (unsigned i = 0; i < nstates-1; i++) {
    for (unsigned j = i+1; j < nstates; j++) {

      unsigned diffs = 0;
      if (_codons[i][0] != _codons[j][0]) {

        diffs++;

      }
      if (_codons[i][1] != _codons[j][1]) {

        diffs++;

      }
      if (_codons[i][2] != _codons[j][2]) {

        diffs++;

      }
      if (diffs == 1) {
        bool synonymous = _amino_acids[i] == _amino_acids[j];
        _Q(i,j) = (synonymous ? 1.0 : omega)*pi[j];
        _Q(j,i) = (synonymous ? 1.0 : omega)*pi[i];
        _Q(i,i) -= _Q(i,j);
        _Q(j,j) -= _Q(j,i);
      }
    }
  }

  double average_rate = 0.0;
  for (unsigned i = 0; i < nstates; i++)
    average_rate -= pi[i]*_Q(i,i);
  double scaling_factor = 3.0/average_rate;
  _Q *= scaling_factor;

  // S is a symmetric matrix
  eigenMatrixXd_t S = eigenMatrixXd_t(_sqrtPi*_Q*_sqrtPiInv);

  // Can use efficient eigensystem solver because S is symmetric
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
  if (solver.info() != Eigen::Success) {

    throw XStrom("Error in the calculation of eigenvectors and eigenvalues of the codon model rate matrix");

  }

  _eigenvectors           = _sqrtPiInv*solver.eigenvectors();
  _inverse_eigenvectors   = solver.eigenvectors().transpose()*_sqrtPi;
  _eigenvalues            = solver.eigenvalues();

}