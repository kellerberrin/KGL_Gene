//
// Created by kellerberrin on 13/12/19.
//

#include "kpl_asrv.h"

#include "kel_distribution.h"

#include <boost/math/distributions/gamma.hpp>

namespace kel = kellerberrin;


namespace kpl = kellerberrin::phylogenetic;


void kpl::ASRV::clear() {

  // Rate homogeneity is the default
  _invar_model = false;
  _ratevar_fixed = false;
  _pinvar_fixed = false;
  _ratevar_ptr = std::make_shared<double>(1.0);
  _pinvar_ptr = std::make_shared<double>(0.0);
  _num_categ = 1;
  recalcASRV();

}





void kpl::ASRV::recalcASRV() {
  // This implementation assumes discrete gamma among-site rate heterogeneity
  // using a _num_categ category discrete gamma distribution with equal category
  // probabilities and Gamma density with mean 1.0 and variance _rate_var.
  // If _invar_model is true, then rate probs will sum to 1 - _pinvar rather than 1
  // and the mean rate will be 1/(1 - _pinvar) rather than 1; the rest of the invariable
  // sites component of the model is handled outside the ASRV class.

  // _num_categ, _rate_var, and _pinvar must all have been assigned in order to compute rates and probs
  if ((not _ratevar_ptr) || (_num_categ == 0) || (not _pinvar_ptr) ) {

    return;

  }

  double pinvar = *_pinvar_ptr;
  assert(pinvar >= 0.0);
  assert(pinvar <  1.0);

  assert(_num_categ > 0);

  double equal_prob = 1.0 /_num_categ;
  double mean_rate_variable_sites = 1.0;

  if (_invar_model) {

    mean_rate_variable_sites = 1.0 / (1.0 - pinvar);

  }

  _rates.assign(_num_categ, mean_rate_variable_sites);
  _probs.assign(_num_categ, equal_prob);

  double rate_variance = *_ratevar_ptr;
  assert(rate_variance >= 0.0);

  if (_num_categ == 1 || rate_variance == 0.0)
    return;

  double alpha = 1.0/rate_variance;
  double beta = rate_variance;

// #define USE_BOOST  1 // Check the timing between the two gamma distribution functions.
#ifdef USE_BOOST
  boost::math::gamma_distribution<> my_gamma(alpha, beta);
  boost::math::gamma_distribution<> my_gamma_plus(alpha + 1.0, beta);
#else
  kel::GammaDistribution k_gamma(alpha, beta);
  kel::GammaDistribution k_gamma_plus(alpha + 1.0, beta);
#endif

  double cum_upper        = 0.0;
  double cum_upper_plus   = 0.0;
  double upper            = 0.0;
  double cum_prob         = 0.0;

  for (unsigned i = 1; i <= _num_categ; ++i) {

    double cum_lower_plus       = cum_upper_plus;
    double cum_lower            = cum_upper;
    cum_prob                    += equal_prob;

    if (i < _num_categ) {

#ifdef USE_BOOST
      upper                   = boost::math::quantile(my_gamma, cum_prob);
      cum_upper_plus          = boost::math::cdf(my_gamma_plus, upper);
      cum_upper               = boost::math::cdf(my_gamma, upper);
#else
      upper                   = k_gamma.quantile(cum_prob);
      cum_upper_plus          = k_gamma_plus.cdf(upper);
      cum_upper               = k_gamma.cdf(upper);
#endif
    }
    else {
      cum_upper_plus          = 1.0;
      cum_upper               = 1.0;
    }

    double numer                = cum_upper_plus - cum_lower_plus;
    double denom                = cum_upper - cum_lower;
    double r_mean               = (denom > 0.0 ? (alpha*beta*numer/denom) : 0.0);
    _rates[i-1]        = r_mean*mean_rate_variable_sites;
  }

}