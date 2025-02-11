//
// Created by kellerberrin on 20/12/19.
//

#include "kpl_mcmc_polytomyprior.h"

#include <limits>
#include <cmath>


namespace kpl = kellerberrin::phylogenetic;


kpl::PolytomyTopoPriorCalculator::PolytomyTopoPriorCalculator() {

  _topo_priors_dirty            = true;
  _is_rooted                    = false;
  _is_resolution_class_prior    = true;
  _C                            = 1.0;
  _ntax                        = 4;
  _counts_dirty                = true;
  _log_scaling_factor          = 10.0;

}

kpl::PolytomyTopoPriorCalculator::~PolytomyTopoPriorCalculator() {

  _counts.clear();
  _topology_prior.clear();
  _nfactors.clear();

}


double kpl::PolytomyTopoPriorCalculator::getLogNormalizedTopologyPrior(unsigned m) {

  if (_topo_priors_dirty) {

    reset();

  }

  assert(m < _topology_prior.size());
  return (_topology_prior[m] - _topology_prior[0]);

}


void kpl::PolytomyTopoPriorCalculator::recalcPriorsImpl() {

  _topology_prior.clear();
  _topology_prior.push_back(0.0);    // This will hold the normalizing constant in the end

  // Figure out the maximum possible value for m, the number of internal nodes
  unsigned maxm = _ntax - (_is_rooted ? 1 : 2);

  if (_is_resolution_class_prior) {
    // _counts vector should have length equal to maxm if everything is ok
    assert(maxm == (unsigned)_counts.size());

    double logC = std::log(_C);
    double log_normalization_constant = 0.0;
    if (_C == 1.0) {

      log_normalization_constant = std::log(maxm);

    }
    else if (_C > 1.0) {
      // factor out largest value to avoid overflow
      double term1  = logC*maxm;
      double term2a = std::exp(-logC*maxm);
      double term2b = std::log(1.0 - term2a);
      double term3  = std::log(_C - 1.0);
      log_normalization_constant = term1 + term2b + term3;

    }
    else {

      log_normalization_constant = std::log(1.0 - std::exp(logC*maxm)) - std::log(1.0 - _C);

    }

    _topology_prior[0] = log_normalization_constant;

    for (unsigned m = 1; m <= maxm; ++m) {

      double logCterm = (double)(maxm - m)*logC;
      double log_count_m = std::log(_counts[m - 1]) + _log_scaling_factor*(double)_nfactors[m - 1];
      double log_v = logCterm - log_count_m;
      _topology_prior.push_back(log_v);

    }

  }
  else {

    double total = 0.0;
    double logC = std::log(_C);
    for (unsigned m = 1; m <= maxm; ++m) {

      double logCterm = (double)(maxm - m)*logC;
      total += std::exp(logCterm);
      _topology_prior.push_back(logCterm);

    }

    _topology_prior[0] = std::log(total);

  }

}

void kpl::PolytomyTopoPriorCalculator::recalcCountsAndPriorsImpl(unsigned n) {

  if (_is_resolution_class_prior) {

    _counts_dirty = true;

  }

  if (_counts_dirty) {

    double scaling_factor = exp(_log_scaling_factor);

    _counts.clear();
    _counts.push_back(1.0); // _counts are always 1 for m = 1

    int last_factor = 0;
    _nfactors.clear();
    _nfactors.push_back(last_factor); // never need to scale this one

    // temporary variables
    double a, b, c;
    const double max_log_double = std::log(std::numeric_limits<double>::max());

    // The value of epsilon is arbitrary, but must be larger than
    // zero and less than scaling_factor
    double epsilon = scaling_factor/10.0;

    // Compute the vector of _counts for the number of internal nodes specified
    // This is the main loop over columns. z is the number of taxa minus 1
    // for rooted trees, and is the number of taxa minus 2 for unrooted trees
    for (unsigned z = 2; z <= n; ++z) {
      // This column is one element longer than the column to its left
      _counts.push_back(0.0);
      _nfactors.push_back(last_factor);

      // _counts[0] is always 1.0 because there is only one star tree topology
      b = _counts[0];

      // This is the loop over rows within the current column.
      // m + 1 is the number of internal nodes.
      for (unsigned m = 1; m < z; ++m) {

        unsigned num_internal_nodes = m + 1;
        double diff = (double)(_nfactors[m - 1] - _nfactors[m]);
        double log_factor = diff*_log_scaling_factor;
        assert(log_factor < max_log_double);
        a = b*exp(log_factor);
        b = _counts[m];
        c = a*((double)(z + num_internal_nodes - 1));

        if (num_internal_nodes < z) {

          c += b*(double)num_internal_nodes;

        }

        if (c > scaling_factor) {

          double incr = floor(log(c - epsilon)/_log_scaling_factor);
          _nfactors[m] += (unsigned)incr;
          last_factor = _nfactors[m];
          _counts[m] = exp(log(c) - incr*_log_scaling_factor);
          b = exp(log(b) - incr*_log_scaling_factor);

        }
        else {

          _counts[m] = c;

        }

      }

    }

    // Now compute the log of the total number of tree topologies
    // over all possible resolution classes (i.e. number of internal nodes)
    // Begin by creating a vector of log _counts and finding the
    // largest value (this will be factored out to avoid overflow)
    std::vector<double> v;
    unsigned sz = (unsigned)_nfactors.size();
    assert(sz == _counts.size());
    double max_log_count = 0.0;

    for (unsigned i = 0; i < sz; ++i) {

      double num_factors = (double)_nfactors[i];
      double log_count = num_factors*_log_scaling_factor + log(_counts[i]);

      if (log_count > max_log_count) {

        max_log_count = log_count;

      }

      v.push_back(log_count);

    }

    // Compute log sum of _counts by factoring out the largest count.
    // Underflow will occur, but only for _counts that are so much
    // smaller than the dominant _counts that the underflow can be
    // ignored for all practical purposes
    double sum = 0.0;
    for (std::vector<double>::const_iterator it = v.begin(); it != v.end(); ++it) {

      double diff = (*it) - max_log_count;
      sum += exp(diff);

    }

    assert(sum > 0.0);

    _log_total_count = log(sum) + max_log_count;
    _counts_dirty = false;

  }
  else {

    _nfactors.clear();
    _counts.clear();

    // _counts_dirty ensures that _counts will be calculated if requested,
    // for example by calling getCount
    _counts_dirty = true;
  }

  // Recalculate the _topology_prior vector too
  recalcPriorsImpl();

  _topo_priors_dirty = false;

}


void kpl::PolytomyTopoPriorCalculator::reset() {

  unsigned num_internal_nodes = (_is_rooted ? (_ntax - 1) : (_ntax - 2));
  recalcCountsAndPriorsImpl(num_internal_nodes);

}

double kpl::PolytomyTopoPriorCalculator::getLogCount(unsigned n, unsigned m) {

  assert((_is_rooted && (m < n)) || (!_is_rooted && (m < n - 1)));

  if (n != _ntax) {

    setNTax(n);

  }

  if (_counts_dirty) {

    reset();

  }

  double nf = (double)(_nfactors[m - 1]);
  double log_count = nf*_log_scaling_factor + log(_counts[m - 1]);
  return log_count;

}

double kpl::PolytomyTopoPriorCalculator::getLogSaturatedCount(unsigned n) {

  if (n != _ntax) {

    setNTax(n);

  }

  if (_counts_dirty) {

    reset();

  }

  unsigned last = (unsigned)(_counts.size() - 1);
  double nf = (double)(_nfactors[last]);
  double log_count = nf*_log_scaling_factor + log(_counts[last]);
  return log_count;

}


std::vector<double> kpl::PolytomyTopoPriorCalculator::getLogCounts() {

  if (_is_resolution_class_prior) {

    _counts_dirty = true;

  }

  if (_counts_dirty) {

    reset();

  }

  std::vector<double> v;
  unsigned sz = _ntax - (_is_rooted ? 0 : 1);
  v.reserve(sz);
  v.push_back(_log_total_count);

  for (unsigned i = 1; i < sz; ++i) {

    double log_Tnm = log(_counts[i - 1]) + _log_scaling_factor*(double)(_nfactors[i - 1]);
    v.push_back(log_Tnm);

  }

  return v;

}


double kpl::PolytomyTopoPriorCalculator::getLogTotalCount(unsigned n) {

  if (n != _ntax) {

    setNTax(n);

  }

  if (_counts_dirty) {

    reset();

  }

  return _log_total_count;

}


std::vector<double> kpl::PolytomyTopoPriorCalculator::getRealizedResClassPriorsVect() {

  if (!_is_resolution_class_prior) {

    _counts_dirty = true;

  }

  if (_topo_priors_dirty || _counts_dirty) {

    reset();

  }

  std::vector<double> v;
  v.reserve(_topology_prior.size());
  v.push_back(0.0);

  unsigned sz = (unsigned)_topology_prior.size();

  // First loop will be to determine largest value, which will be factored out
  // the second time through so that the total does not overflow
  double log_factored_out = 0.0;
  for (unsigned i = 1; i < sz; ++i) {

    double c = _counts[i - 1];
    double nf = (double)_nfactors[i - 1];
    double log_Tnm = log(c) + _log_scaling_factor*nf;
    double log_prior = log_Tnm + _topology_prior[i];

    v.push_back(log_prior);

    if (log_prior > log_factored_out) {

      log_factored_out = log_prior;

    }

  }

  // Now we can compute the total
  double total = 0.0;
  std::vector<double>::const_iterator it = v.begin();

  for (++it; it != v.end(); ++it) {

    total += exp((*it) - log_factored_out);

  }

  v[0] = log(total) + log_factored_out;

  return v;

}


unsigned kpl::PolytomyTopoPriorCalculator::sample(Lot::SharedPtr rng) {

  std::vector<double> v = getRealizedResClassPriorsVect();
  double u = rng->uniform();
  double z = v[0];
  double cum = 0.0;

  for (unsigned i = 1; i < v.size(); ++i) {

    cum += exp(v[i] - z);

    if (u <= cum) {

      return i;

    }

  }

  assert(0);

  return (unsigned)(v.size() - 1);

}


void kpl::PolytomyTopoPriorCalculator::setNTax(unsigned new_ntax) {

  if (_ntax != new_ntax) {

    // Set _ntax to the new value
    assert(new_ntax > (unsigned)(_is_rooted ? 1 : 2));
    _ntax = new_ntax;

    _counts_dirty = true;
    _topo_priors_dirty = true;

  }

}


unsigned kpl::PolytomyTopoPriorCalculator::getNTax() const {

  return _ntax;

}


void kpl::PolytomyTopoPriorCalculator::chooseRooted() {

  if (!_is_rooted) {

    _is_rooted = true;
    _topo_priors_dirty = true;

  }

}


void kpl::PolytomyTopoPriorCalculator::chooseUnrooted() {

  if (_is_rooted) {

    _is_rooted = false;
    _topo_priors_dirty = true;

  }

}


void kpl::PolytomyTopoPriorCalculator::chooseResolutionClassPrior() {

  if (!_is_resolution_class_prior) {

    _is_resolution_class_prior = true;
    _topo_priors_dirty = true;

  }

}


void kpl::PolytomyTopoPriorCalculator::choosePolytomyPrior() {

  if (_is_resolution_class_prior) {

    _is_resolution_class_prior = false;
    _topo_priors_dirty = true;

  }

}


void kpl::PolytomyTopoPriorCalculator::setC(double c) {

  assert(c > 0.0);
  if (c != _C) {

    _C = c;
    _topo_priors_dirty = true;

  }

}

double kpl::PolytomyTopoPriorCalculator::getC() const {

  return _C;

}


void kpl::PolytomyTopoPriorCalculator::setLogScalingFactor(double lnf) {

  assert(lnf > 0.0);

  if (lnf != _log_scaling_factor) {

    _log_scaling_factor = lnf;
    _counts_dirty = true;
    _topo_priors_dirty = true;

  }

}


double kpl::PolytomyTopoPriorCalculator::getLogScalingFactor() const {

  return _log_scaling_factor;

}


std::vector<double> kpl::PolytomyTopoPriorCalculator::getCountsVect() {

  if (_counts_dirty) {

    reset();

  }

  return _counts;

}


std::vector<int> kpl::PolytomyTopoPriorCalculator::getNFactorsVect() {

  if (_counts_dirty) {

    reset();

  }

  return _nfactors;

}


std::vector<double> kpl::PolytomyTopoPriorCalculator::getTopoPriorVect() {
  //@POL this function could be const were it not for lazy evaluation
  if (_topo_priors_dirty) {

    reset();

  }
  return _topology_prior;

}


double kpl::PolytomyTopoPriorCalculator::getLogTopoProb(std::shared_ptr<Tree> t) {

  unsigned n = t->numLeaves();
  assert(n > 3);

  if (n != getNTax()) {

    setNTax(n);

  }

  unsigned m = t->numInternals();
  double topo_prior = getLogTopologyPrior(m);
  return topo_prior;

}


double kpl::PolytomyTopoPriorCalculator::getLogTopologyPrior(unsigned m) {

  if (_topo_priors_dirty) {

    reset();

  }

  assert(m < _topology_prior.size());

  return _topology_prior[m];

}


double kpl::PolytomyTopoPriorCalculator::getLogNormConstant() {

  if (_topo_priors_dirty) {

    reset();

  }

  return _topology_prior[0];

}


bool kpl::PolytomyTopoPriorCalculator::isResolutionClassPrior() const {

  return _is_resolution_class_prior;

}


bool kpl::PolytomyTopoPriorCalculator::isPolytomyPrior() const {

  return !_is_resolution_class_prior;

}


bool kpl::PolytomyTopoPriorCalculator::isRooted() const {

  return _is_rooted;

}


bool kpl::PolytomyTopoPriorCalculator::isUnrooted() const {

  return !_is_rooted;

}