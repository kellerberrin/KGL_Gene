//
// Created by kellerberrin on 13/05/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_deconvolv_app.h"
#include "kgd_ibdpath.h"



namespace kgd = kellerberrin::deconvolv;



void kgd::IBDpath::init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator) {

  ibd_random_generator_ = randomGenerator;
  setNLoci(dEploidIO.nLoci());
  setKstrain(dEploidIO.kStrain());
  setTheta(1.0 / (double) kStrain());

  IBD_path_change_at_ = std::vector<double>(nLoci());

  // compute likelihood surface
  ibd_prob_cache_ = std::make_shared<const SiteProbabilityCache>(dEploidIO.getAltCount(),
                                                                 dEploidIO.getRefCount(),
                                                                 dEploidIO.ibdParameters().betaBinomialConstant(), /*  double scalingConst */
                                                                 dEploidIO.ibdParameters().baseCountError(), /* double err */
                                                                 dEploidIO.ibdParameters().cacheGridSize()); /* cache size */

  // initialize haplotype prior
  h_prior_.initializeHprior(kStrain(), dEploidIO.getPlaf());

  makeIbdTransProbs();

  // initialize fm_
  f_sum_state_ = std::vector<double>(h_prior_.nStates());

  // initialize ibd_configure_path_
  ibd_configure_path_ = std::vector<size_t>(nLoci());

  // initialize recombination probabilities;
  ibd_recomb_probs_ = std::make_shared<const IBDRecombProbs>(dEploidIO.getPosition(),
                                                             dEploidIO.nLoci(),
                                                             dEploidIO.ibdParameters().averageCentimorganDistance(),
                                                             dEploidIO.ibdParameters().parameterG(),
                                                             dEploidIO.useConstRecomb(),
                                                             dEploidIO.ibdParameters().constRecombProb());

  current_IBD_path_change_at_ = std::vector<double>(nLoci());


  computeUniqueEffectiveKCount();

}


void kgd::IBDpath::McmcUpdateStep(const std::vector<double>& CurrentProportion) {

  std::vector<double> effectiveKPrior = computeEffectiveKPrior(theta());

  std::vector<double> statePrior = computeStatePrior(effectiveKPrior);

  ExecEnv::log().info("effectiveKPrior.size: {}, statePrior.size: {}, theta: {}", effectiveKPrior.size(), statePrior.size(), theta());

  /// First build the path likelihood
  computeIbdPathFwdProb(CurrentProportion, statePrior);
  /// Now sample path given matrix
  ibdSamplePath(statePrior);
  /// Compute new theta after all proportion and haplotypes are up to date.
  computeAndUpdateTheta();

}


double kgd::IBDpath::UpdateHaplotypesFromPrior(size_t strain, size_t loci) const {

  size_t config_idx = idbConfigurePath()[loci];

  auto set_value = getHprior().gethSet()[config_idx][strain];

  return static_cast<double>(set_value);

}


void kgd::IBDpath::ibdSamplePath(const std::vector<double>& statePrior) {

  int lociIdx = nLoci() - 1;

  std::vector<double> tmpProp = fm_[lociIdx];

  Utility::normalizeBySum(tmpProp);

  ibd_configure_path_[lociIdx] = Utility::sampleIndexGivenProp(ibd_random_generator_, tmpProp);

  assert(fm_.size() == nLoci());

  while (lociIdx > 0) {

    lociIdx--;

    std::vector<double> vNoRecomb = Utility::vecProd(ibd_trans_probs_[h_prior_.getStateIdx()[ibd_configure_path_[lociIdx + 1]]], fm_[lociIdx]);

    assert(vNoRecomb.size() == h_prior_.nStateEntries());

    std::vector<double> vRecomb = fm_[lociIdx];

    assert(vRecomb.size() == h_prior_.nStateEntries());

    std::vector<double> prop(h_prior_.nStateEntries());

    for (size_t i = 0; i < prop.size(); i++) {

      prop[i] = vNoRecomb[i] * ibd_recomb_probs_->getNoRecombProbAtLoci(lociIdx) +
                vRecomb[i] * ibd_recomb_probs_->getRecombProbAtLoci(lociIdx) * statePrior[ibd_configure_path_[lociIdx + 1]];

    }

    tmpProp = prop;
    Utility::normalizeBySum(tmpProp);
    ibd_configure_path_[lociIdx] = Utility::sampleIndexGivenProp(ibd_random_generator_, tmpProp);

    assert(ibd_configure_path_[lociIdx] < h_prior_.nStateEntries());

  }

}


std::vector<size_t> kgd::IBDpath::findAllIndex(const std::vector<size_t> &index_array, size_t index) {

  std::vector<size_t> ret;

  for (size_t i = 0; i < index_array.size(); i++) {

    if (index_array[i] == index) {

      ret.push_back(i);

    }

  }

  return ret;

}


void kgd::IBDpath::buildPathProbabilityForPainting(const std::vector<double>& proportion) {

  //vector <double> effectiveKPrior = computeEffectiveKPrior(theta());
  std::vector<double> effectiveKPrior = std::vector<double>(h_prior_.nStates(), 1.0 / h_prior_.nStates());
  std::vector<double> statePrior = computeStatePrior(effectiveKPrior);

  // First building the path likelihood
  computeIbdPathFwdProb(proportion, statePrior);

  // Reshape Fwd
  std::vector<std::vector<double>> reshapedFwd = reshapeProbs(fm_);

  computeIbdPathBwdProb(proportion, effectiveKPrior, statePrior);
  // Reshape Bwd
  std::vector<std::vector<double>> reshapedBwd = reshapeProbs(bwd_);

  // Combine Fwd Bwd
  combineFwdBwd(reshapedFwd, reshapedBwd);

}


void kgd::IBDpath::computeIbdPathBwdProb(const std::vector<double>& proportion,
                                         const std::vector<double>& effectiveKPrior,
                                         const std::vector<double>& statePrior) {

  //# assuming each ibd state has equal probabilities, transform it into ibd configurations
  //dout << " start building ibd bwd_ "<< endl;
  std::vector<double> tmp = std::vector<double>(h_prior_.getStateIdxFreq().size());

  assert(effectiveKPrior.size() == h_prior_.getStateIdxFreq().size());

  for (size_t i = 0; i < tmp.size(); i++) {

    tmp[i] = effectiveKPrior[i] / (double) h_prior_.getStateIdxFreq()[i];

  }

  std::vector<double> tmpBw = std::vector<double>(h_prior_.nStateEntries());

  for (size_t j = 0; j < tmpBw.size(); j++) {

    for (size_t i = 0; i < tmp.size(); i++) {

      tmpBw[j] += tmp[i] * ibd_trans_probs_[i][j];

    }

  }

  bwd_.push_back(tmpBw);

  for (size_t rev_site_i = 1; rev_site_i < nLoci(); ++rev_site_i) {

    size_t site_i = nLoci() - rev_site_i;

    std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, site_i);

    //vector<double> lk = vector <double> (h_prior_.nState(), 1.0);
    std::vector<double> bSumState = std::vector<double>(h_prior_.nStates());
    for (size_t i = 0; i < bSumState.size(); i++) {

      for (size_t j = 0; j < h_prior_.nStateEntries(); j++) {

        bSumState[i] += ibd_trans_probs_[i][j] * bwd_.back()[j];

      }

    }

    std::vector<double> vNoRecomb(h_prior_.nStateEntries());
    for (size_t i = 0; i < h_prior_.getStateIdx().size(); i++) {

      vNoRecomb[i] = bSumState[h_prior_.getStateIdx()[i]];

    }

    for (size_t state = 0; state < h_prior_.nStateEntries(); ++state) {

      tmpBw[state] = 0;

      for (size_t j = 0; j < lk.size(); j++) {

        tmpBw[state] += (lk[j] * bwd_.back()[j]) * ibd_recomb_probs_->getRecombProbAtLoci(site_i - 1);

      }

      tmpBw[state] *= statePrior[state];
      tmpBw[state] += lk[state] * (ibd_recomb_probs_->getNoRecombProbAtLoci(site_i - 1)) * vNoRecomb[state];
      tmpBw[state] *= h_prior_.getPriorProb()[site_i][state];

    }

    Utility::normalizeBySum(tmpBw);

    bwd_.push_back(tmpBw);

  }

  std::reverse(bwd_.begin(), bwd_.end());

}


void kgd::IBDpath::computeIbdPathFwdProb(const std::vector<double>& proportion, const std::vector<double>& statePrior) {

  fm_.clear();

  std::vector<double> vPrior = Utility::vecProd(statePrior, h_prior_.getPriorProb()[0]);

  std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, 0);

  updateFmAtSiteI(vPrior, lk);

  for (size_t site_i = 1; site_i < nLoci(); ++site_i) {

    std::vector<double> vNoRec;

    for (size_t stateIdxTmp : h_prior_.getStateIdx()) {

      vNoRec.push_back(f_sum_state_[stateIdxTmp]);

    }

    for (size_t state = 0; state < h_prior_.nStateEntries(); state++) {

      vPrior[state] = ((vNoRec[state] * ibd_recomb_probs_->getNoRecombProbAtLoci(site_i))
                   + (f_sum_ * ibd_recomb_probs_->getRecombProbAtLoci(site_i) * statePrior[state]))
                  * h_prior_.getPriorProb()[site_i][state];

    }

    lk = computeLlkOfStatesAtSiteI(proportion, site_i);

    updateFmAtSiteI(vPrior, lk);

  }

}


void kgd::IBDpath::updateFmAtSiteI(const std::vector<double> &prior, const std::vector<double> &llk) {

  std::vector<double> postAtSiteI = Utility::vecProd(prior, llk);
  //normalizeByMax(postAtSiteI);

  Utility::normalizeBySum(postAtSiteI);

  fm_.push_back(postAtSiteI);

  f_sum_ = Utility::sumOfVec(postAtSiteI);

  for (size_t i = 0; i < f_sum_state_.size(); i++) {

    f_sum_state_[i] = 0;

    for (size_t j = 0; j < h_prior_.nStateEntries(); j++) {

      f_sum_state_[i] += ibd_trans_probs_[i][j] * postAtSiteI[j];

    }

  }

}


double kgd::IBDpath::bestPath(const std::vector<double>& proportion, double err) const {

  double sumLLK = 0.0;

  for (size_t siteI = 0; siteI < nLoci(); siteI++) {

    std::vector<double> tmp;

    for (size_t j = 0; j < fm_[siteI].size(); j++) {

      tmp.push_back(std::exp(std::log(fm_[siteI][j]) + std::log(bwd_[siteI][j])));

    }

    Utility::normalizeBySum(tmp);

    size_t indx = std::distance(tmp.begin(), std::max_element(tmp.begin(), tmp.end()));

    std::vector<size_t> hSetI = h_prior_.gethSet()[indx];

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += static_cast<double>(hSetI[j]) * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    if ((qs > 0) & (qs < 1)) {

      sumLLK += siteLogBetaLLK(siteI, qs2);

    }

  }

  return sumLLK;

}


void kgd::IBDpath::combineFwdBwd(const std::vector<std::vector<double>> &reshapedFwd, const std::vector<std::vector<double>> &reshapedBwd) {

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> tmp;

    for (size_t j = 0; j < reshapedFwd[i].size(); j++) {

      tmp.push_back(exp(log(reshapedFwd[i][j]) + log(reshapedBwd[i][j])));

    }

    Utility::normalizeBySum(tmp);

    fwdbwd_.push_back(tmp);

  }

}

// This creates a matrix rows = states (52), columns = state_entries (454).
// For each state (row) only the entries that correspond to that state are set to 1.0, else 0.0.
void kgd::IBDpath::makeIbdTransProbs() {

  assert(ibd_trans_probs_.size() == 0);

  for (size_t i = 0; i < h_prior_.nStates(); i++) {

    std::vector<double> transProbRow(h_prior_.nStateEntries());
    std::vector<size_t> wi = findAllIndex(h_prior_.getStateIdx(), i);

    ExecEnv::log().info("IBDPath; transProbRow.size: {}, wi.size: {}", transProbRow.size(), wi.size());

    for (size_t wii : wi) {

      transProbRow[wii] = 1;

    }

    ibd_trans_probs_.push_back(transProbRow);

  }

}


std::vector<std::string> kgd::IBDpath::getIBDprobsHeader() const {

  return h_prior_.getIBDconfigureHeader();

}


std::vector<std::vector<double> > kgd::IBDpath::reshapeProbs(const std::vector<std::vector<double> > &probs) const {

  assert(nLoci() == probs.size());

  std::vector<std::vector<double> > ret;

  for (size_t siteIndex = 0; siteIndex < nLoci(); siteIndex++) {

    size_t previousStateIdx = 0;
    std::vector<double> tmpRow;
    double cumProb = 0;

    for (size_t prob_ij = 0; prob_ij < probs[siteIndex].size(); prob_ij++) {

      cumProb += probs[siteIndex][prob_ij];

      if (previousStateIdx != h_prior_.getStateIdx()[prob_ij]) {
        cumProb -= probs[siteIndex][prob_ij];
        previousStateIdx++;
        tmpRow.push_back(cumProb);
        cumProb = probs[siteIndex][prob_ij];
      }
    }

    tmpRow.push_back(cumProb);
    Utility::normalizeBySum(tmpRow);
    ret.push_back(tmpRow);

  }

  return ret;

}


std::vector<double> kgd::IBDpath::computeEffectiveKPrior(double theta) {
  //#Calculate state prior given theta (theta is prob IBD)

  std::vector<double> binomial_vector;

  for (size_t i = 0; i < kStrain(); ++i) {

    double binomial_i = Utility::binomialPdf(i, (kStrain() - 1), theta);

    binomial_vector.push_back(binomial_i);

  }

  std::vector<double> effectiveKPrior;

  for (auto effectiveKtmp : h_prior_.getEffectiveK()) {

    size_t effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx < kStrain());

    double  binomial_ratio = binomial_vector[effectiveKidx] / static_cast<double>(unique_effectiveK_count_[effectiveKidx]);

    effectiveKPrior.push_back(binomial_ratio);

  }

  return effectiveKPrior;

}


std::vector<double> kgd::IBDpath::computeStatePrior(const std::vector<double>& effectiveKPrior) {

  std::vector<double> ret;

  for (size_t stateIdxTmp : h_prior_.getStateIdx()) {

    ret.push_back(effectiveKPrior[stateIdxTmp]);

  }

  return ret;

}


void kgd::IBDpath::computeAndUpdateTheta() {

  std::vector<size_t> obsState;

  size_t previousState = 0;
  size_t atSiteI = 0;

  for (size_t a : ibd_configure_path_) {

    if (a != previousState) {

      obsState.push_back(a);

    }

    if (h_prior_.getStateIdx()[a] != h_prior_.getStateIdx()[previousState]) {

      IBD_path_change_at_[atSiteI] += 1.0;
      current_IBD_path_change_at_[atSiteI] = 1.0;

    } else {

      current_IBD_path_change_at_[atSiteI] = 0.0;

    }

    previousState = a;
    atSiteI++;

  }

  size_t sumOfKeffStates = 0;
  size_t sccs = 0;

  for (size_t obs : obsState) {

    sumOfKeffStates += h_prior_.getEffectiveK()[obs] - 1;
    sccs += kStrain() - h_prior_.getEffectiveK()[obs];

  }

  setTheta(Utility::rBeta(sccs + 1.0, sumOfKeffStates + 1.0, ibd_random_generator_));

}


void kgd::IBDpath::computeUniqueEffectiveKCount() {

  unique_effectiveK_count_ = std::vector<size_t>(kStrain());

  for (size_t effectiveKtmp : h_prior_.getEffectiveK()) {

    int effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx >= 0);

    unique_effectiveK_count_[effectiveKidx]++;

  }

}


std::vector<double> kgd::IBDpath::computeLlkOfStatesAtSiteI(const std::vector<double>& proportion, size_t siteI, double err) {

  std::vector<double> llks;

  for (std::vector<size_t> hSetI : h_prior_.gethSet()) {

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += static_cast<double>(hSetI[j]) * proportion[j];

    }

    double qs2 = qs * (1 - err) + (1 - qs) * err;

    double logBetaLLK = siteLogBetaLLK(siteI, qs2);

    llks.push_back(logBetaLLK);

  }

  double maxllk = Utility::max_value(llks);

  std::vector<double> ret;

  for (double llk : llks) {

    double normalized = std::exp(llk - maxllk);

    if (normalized == 0) {

      normalized = std::numeric_limits<double>::min();

    }

    ret.push_back(normalized);

  }

  return ret;
}


