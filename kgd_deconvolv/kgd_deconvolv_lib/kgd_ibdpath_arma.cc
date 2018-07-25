//
// Created by kellerberrin on 19/07/18.
//


#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "kgd_deconvolv_app.h"
#include "kgd_ibdpath_arma.h"



namespace kgd = kellerberrin::deconvolv;



void kgd::IBDPath::init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator) {

  ibd_random_generator_ = randomGenerator;
  setNLoci(dEploidIO.nLoci());
  setKstrain(dEploidIO.kStrain());
  setTheta(1.0 / (double) kStrain());

  IBD_path_change_at_ = std::vector<double>(nLoci());

  // compute likelihood surface
  ibd_prob_cache_ = std::make_shared<const SiteProbabilityCache>(dEploidIO.getAltCount(),
                                                                 dEploidIO.getRefCount(),
                                                                 100.0, /*  double scalingConst */
                                                                 0.01, /* double err */
                                                                 999 /* cache size */ );

  // initialize haplotype prior
  h_prior_.initializeHprior(kStrain(), dEploidIO.getPlaf());

  makeIbdTransProbs();

  // initialize fm_
  f_sum_state_ = std::vector<double>(h_prior_.nPattern());

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


void kgd::IBDPath::McmcUpdateStep(const std::vector<double>& CurrentProportion) {

  arma::rowvec effectiveKPrior;
  computeEffectiveKPrior(effectiveKPrior, theta());

  arma::rowvec statePrior;
  computeStatePrior(statePrior, effectiveKPrior);
  /// First build the path likelihood
  computeIbdPathFwdProb(CurrentProportion, statePrior);
  /// Now sample path given matrix
  ibdSamplePath(statePrior);
  /// Compute new theta after all proportion and haplotypes are up to date.
  computeAndUpdateTheta();

}


double kgd::IBDPath::UpdateHaplotypesFromPrior(size_t strain, size_t loci) const {

  size_t config_idx = idbConfigurePath()[loci];

  auto set_value = getHprior().gethSet()[config_idx][strain];

  return static_cast<double>(set_value);

}


void kgd::IBDPath::ibdSamplePath(const arma::rowvec& statePrior) {

  size_t lociIdx = nLoci() - 1;

  // Find the max strain value from fm_(loci, strains)
  double max_value = 0.0;
  size_t max_index = 0;
  for (size_t i = 0; i < fm_.n_cols; ++i) {

    if (fm_(lociIdx, i) > max_value) {

      max_value = fm_(lociIdx, i);
      max_index = i;

    }

  }

  ibd_configure_path_(lociIdx) = max_index;

  assert(fm_.n_rows == nLoci());

  while (lociIdx > 0) {

    lociIdx--;

    arma::vec vNoRecomb = ibd_trans_probs_(h_prior_.getStateIdx()[ibd_configure_path_(lociIdx + 1)]) % fm_(lociIdx);

    assert(vNoRecomb.n_rows == h_prior_.nState());

    for (size_t i = 0; i < h_prior_.nState(); i++) {

      prop  vNoRecomb[i] * ibd_recomb_probs_->getNoRecombProbAtLoci(lociIdx) +
                vRecomb[i] * ibd_recomb_probs_->getRecombProbAtLoci(lociIdx) * statePrior[ibd_configure_path_[lociIdx + 1]];

    }

    tmpProp = prop;
    Utility::normalizeBySum(tmpProp);
    ibd_configure_path_[lociIdx] = Utility::sampleIndexGivenProp(ibd_random_generator_, tmpProp);

    assert(ibd_configure_path_[lociIdx] < h_prior_.nState());

  }

}


std::vector<size_t> kgd::IBDPath::findWhichIsSomething(const std::vector<size_t>& tmpOp, size_t something) {

  std::vector<size_t> ret;

  for (size_t i = 0; i < tmpOp.size(); i++) {

    if (tmpOp[i] == something) {
      ret.push_back(i);
    }

  }

  return ret;

}


void kgd::IBDPath::buildPathProbabilityForPainting(const std::vector<double>& proportion) {

  //vector <double> effectiveKPrior = computeEffectiveKPrior(theta());
  std::vector<double> effectiveKPrior = std::vector<double>(h_prior_.nPattern(), 1.0 / h_prior_.nPattern());
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


void kgd::IBDPath::computeIbdPathBwdProb(const std::vector<double>& proportion,
                                         const std::vector<double>& effectiveKPrior,
                                         const std::vector<double>& statePrior) {

  //# assuming each ibd state has equal probabilities, transform it into ibd configurations
  //dout << " start building ibd bwd_ "<< endl;
  std::vector<double> tmp = std::vector<double>(h_prior_.getStateIdxFreq().size());

  assert(effectiveKPrior.size() == h_prior_.getStateIdxFreq().size());

  for (size_t i = 0; i < tmp.size(); i++) {

    tmp[i] = effectiveKPrior[i] / (double) h_prior_.getStateIdxFreq()[i];

  }

  std::vector<double> tmpBw = std::vector<double>(h_prior_.nState());

  for (size_t j = 0; j < tmpBw.size(); j++) {

    for (size_t i = 0; i < tmp.size(); i++) {

      tmpBw[j] += tmp[i] * ibd_trans_probs_[i][j];

    }

  }

  bwd_.push_back(tmpBw);

  for (size_t rev_siteI = 1; rev_siteI < nLoci(); rev_siteI++) {

    size_t siteI = nLoci() - rev_siteI;

    std::vector<double> lk = computeLlkOfStatesAtSiteI(proportion, siteI);

    //vector<double> lk = vector <double> (h_prior_.nState(), 1.0);
    std::vector<double> bSumState = std::vector<double>(h_prior_.nPattern());
    for (size_t i = 0; i < bSumState.size(); i++) {

      for (size_t j = 0; j < h_prior_.nState(); j++) {

        bSumState[i] += ibd_trans_probs_[i][j] * bwd_.back()[j];

      }

    }

    std::vector<double> vNoRecomb(h_prior_.nState());
    for (size_t i = 0; i < h_prior_.getStateIdx().size(); i++) {

      vNoRecomb[i] = bSumState[h_prior_.getStateIdx()[i]];

    }

    for (size_t i = 0; i < h_prior_.nState(); i++) {

      tmpBw[i] = 0;

      for (size_t j = 0; j < lk.size(); j++) {

        tmpBw[i] += (lk[j] * bwd_.back()[j]) * ibd_recomb_probs_->getRecombProbAtLoci(siteI - 1);

      }

      tmpBw[i] *= statePrior[i];
      tmpBw[i] += lk[i] * (ibd_recomb_probs_->getNoRecombProbAtLoci(siteI - 1)) * vNoRecomb[i];
      tmpBw[i] *= h_prior_.getPriorProb()[i][siteI];

    }

    Utility::normalizeBySum(tmpBw);
    bwd_.push_back(tmpBw);

  }

  reverse(bwd_.begin(), bwd_.end());

}


void kgd::IBDPath::computeIbdPathFwdProb(const std::vector<double>& proportion, const arma::rowvec& statePrior) {

  fm_.reset();

  arma::rowvec priorProbTrans(h_prior_.getPriorProbTrans()[0]);

  arma::rowvec lk;
  computeLlkOfStatesAtSiteI(lk, proportion, 0);

  updateFmAtSiteI((statePrior % priorProbTrans), lk);

  for (size_t siteI = 1; siteI < nLoci(); siteI++) {

    std::vector<double> vNoRec;

    for (size_t stateIdxTmp : h_prior_.getStateIdx()) {

      vNoRec.push_back(f_sum_state_[stateIdxTmp]);

    }

    for (size_t i = 0; i < h_prior_.nState(); i++) {

      vPrior[i] = ((vNoRec[i] * ibd_recomb_probs_->getNoRecombProbAtLoci(siteI))
                   + (f_sum_ * ibd_recomb_probs_->getRecombProbAtLoci(siteI) * statePrior[i]))
                  * h_prior_.getPriorProbTrans()[siteI][i];

    }

    lk = computeLlkOfStatesAtSiteI(proportion, siteI);

    updateFmAtSiteI(vPrior, lk);

  }

}


void kgd::IBDPath::updateFmAtSiteI(const arma::rowvec &prior, const arma::rowvec &llk) {

  arma::rowvec postAtSiteI(prior % llk);

  double vecsum = arma::accu(postAtSiteI);

  postAtSiteI = postAtSiteI / vecsum;

  fm_.push_back(postAtSiteI);

  f_sum_ = Utility::sumOfVec(postAtSiteI);

  for (size_t i = 0; i < f_sum_state_.size(); i++) {

    f_sum_state_[i] = 0;

    for (size_t j = 0; j < h_prior_.nState(); j++) {

      f_sum_state_[i] += ibd_trans_probs_[i][j] * postAtSiteI[j];

    }

  }

}


double kgd::IBDPath::bestPath(const std::vector<double>& proportion, double err) const {

  double sumLLK = 0.0;

  for (size_t siteI = 0; siteI < nLoci(); siteI++) {

    std::vector<double> tmp;

    for (size_t j = 0; j < fm_[siteI].size(); j++) {

      tmp.push_back(std::exp(std::log(fm_[siteI][j]) + std::log(bwd_[siteI][j])));

    }

    Utility::normalizeBySum(tmp);

    size_t indx = std::distance(tmp.begin(), std::max_element(tmp.begin(), tmp.end()));

    std::vector<int> hSetI = h_prior_.gethSet()[indx];

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


void kgd::IBDPath::combineFwdBwd(const std::vector<std::vector<double>> &reshapedFwd, const std::vector<std::vector<double>> &reshapedBwd) {

  for (size_t i = 0; i < nLoci(); i++) {

    std::vector<double> tmp;

    for (size_t j = 0; j < reshapedFwd[i].size(); j++) {

      tmp.push_back(exp(log(reshapedFwd[i][j]) + log(reshapedBwd[i][j])));

    }

    Utility::normalizeBySum(tmp);

    fwdbwd_.push_back(tmp);

  }

}


void kgd::IBDPath::makeIbdTransProbs() {

  assert(ibd_trans_probs_.size() == 0);

  for (size_t i = 0; i < h_prior_.nPattern(); i++) {

    std::vector<double> transProbRow(h_prior_.nState());
    std::vector<size_t> wi = findWhichIsSomething(h_prior_.getStateIdx(), i);

    for (size_t wii : wi) {

      transProbRow[wii] = 1;

    }

    ibd_trans_probs_.push_back(transProbRow);

  }

}


std::vector<std::string> kgd::IBDPath::getIBDprobsHeader() const {

  return h_prior_.getIBDconfigureHeader();

}


std::vector<std::vector<double> > kgd::IBDPath::reshapeProbs(const std::vector<std::vector<double> > &probs) const {

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


void kgd::IBDPath::computeEffectiveKPrior(arma::rowvec& effectiveKPrior, double theta) {
  //#Calculate state prior given theta (theta is prob IBD)

  arma::rowvec binomial_vector;

  for (size_t i = 0; i < kStrain(); ++i) {

    double binomial_i = Utility::binomialPdf( i, (kStrain() - 1), theta);
    binomial_vector << binomial_i;

  }

  effectiveKPrior.reset();

  for (auto effectiveKtmp : h_prior_.getEffectiveK()) {

    size_t effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx < kStrain());

    double  binomial_ratio = binomial_vector[effectiveKidx] / static_cast<double>(unique_effectiveK_count_[effectiveKidx]);

    effectiveKPrior << binomial_ratio;

  }

}


void kgd::IBDPath::computeStatePrior(arma::rowvec& statePrior, const arma::rowvec& effectiveKPrior) {

  statePrior.reset();

  for (size_t stateIdxTmp : h_prior_.getStateIdx()) {

    statePrior << effectiveKPrior(stateIdxTmp);

  }

}


void kgd::IBDPath::computeAndUpdateTheta() {

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


void kgd::IBDPath::computeUniqueEffectiveKCount() {

  unique_effectiveK_count_ = std::vector<int>(kStrain());

  for (size_t effectiveKtmp : h_prior_.getEffectiveK()) {

    int effectiveKidx = effectiveKtmp - 1;

    assert(effectiveKidx >= 0);

    unique_effectiveK_count_[effectiveKidx]++;

  }

}


void kgd::IBDPath::computeLlkOfStatesAtSiteI(arma::rowvec& lk, const std::vector<double>& proportion, size_t siteI, double err) {

  lk.reset();

  arma::rowvec llks;

  for (std::vector<int> hSetI : h_prior_.gethSet()) {

    double qs = 0;

    for (size_t j = 0; j < kStrain(); j++) {

      qs += static_cast<double>(hSetI[j]) * proportion[j];

    }

    double qs2 = qs + (1 - (2 * qs)) * err;

    double logBetaLLK = siteLogBetaLLK(siteI, qs2);

    llks << logBetaLLK;

  }

  double maxllk = llks.max();

  for (auto llk : llks) {

    double normalized = std::exp(llk - maxllk);

    if (normalized == 0) {

      normalized = std::numeric_limits<double>::min();

    }

    lk << normalized;

  }

}


