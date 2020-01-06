//
// Created by kellerberrin on 13/12/19.
//

#include "kpl_qmatrix.h"

#include <numeric>

namespace kpl = kellerberrin::phylogenetic;

// member function bodies go here


kpl::QMatrix::QMatrix() {
  //std::cout << "Creating a QMatrix object" << std::endl;
}


kpl::QMatrix::~QMatrix() {
  //std::cout << "Destroying a QMatrix object" << std::endl;
}


void kpl::QMatrix::setActive(bool activate) {

  _is_active = activate;
  recalcRateMatrix();

}


void kpl::QMatrix::clear() {

  _is_active = false;
  _state_freqs_fixed = false;
  _exchangeabilities_fixed = false;
  _omega_fixed = false;

}


void kpl::QMatrix::fixStateFreqs(bool is_fixed) {

  _state_freqs_fixed = is_fixed;

}


void kpl::QMatrix::fixExchangeabilities(bool is_fixed) {

  _exchangeabilities_fixed = is_fixed;

}


void kpl::QMatrix::fixOmega(bool is_fixed) {

  _omega_fixed = is_fixed;

}


bool kpl::QMatrix::isFixedStateFreqs() const {

  return _state_freqs_fixed;

}


bool kpl::QMatrix::isFixedExchangeabilities() const {

  return _exchangeabilities_fixed;

}


bool kpl::QMatrix::isFixedOmega() const {

  return _omega_fixed;

}


void kpl::QMatrix::normalizeFreqsOrExchangeabilities(QMatrix::freq_xchg_ptr_t v) {
  // Be sure elements of v sum to 1.0 and assert that they are all positive
  double sum_v = std::accumulate(v->begin(), v->end(), 0.0);

  for (auto & x : *v) {

    assert(x > 0.0);
    x /= sum_v;

  }

}


