/*
 * kgd_deploid is used for deconvoluting Plasmodium falciparum genome from
 * mix-infected patient sample.
 *
 * Copyright (C) 2016-2017 University of Oxford
 *
 * Author: Sha (Joe) Zhu
 *
 * This file is part of kgd_deploid.
 *
 * kgd_deploid is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */



#ifndef KGD_IBD_H
#define KGD_IBD_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_mersenne_twister.h"
#include "kgd_dEploidIO.h"



namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



int nchoose2(int n);

bool twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2);

std::vector<std::vector<int> > unique(std::vector<std::vector<int> > &mat);

std::vector<int> convertIntToBinary(int x, size_t len);

std::vector<std::vector<int> > enumerateBinaryMatrixOfK(size_t k);

// The IBDconfiguration is used for index, which should be non-negative, use int, any thing below zero should throw.
class IBDconfiguration {
#ifdef UNITTEST
  friend class TestIBDconfig;
  friend class TestHprior;
#endif

  friend class Hprior;

  friend class McmcMachinery;

  IBDconfiguration() = default;
  ~IBDconfiguration() = default;

  void buildIBDconfiguration(size_t k = 5);

  size_t kStrain_;

  void setKstrain(const size_t setTo) { this->kStrain_ = setTo; }

  size_t kStrain() const { return this->kStrain_; }

  std::vector<std::vector<int> > op;
  std::vector<std::vector<int> > pairToEmission;
  std::vector<std::vector<size_t> > pairList;
  std::vector<std::vector<int> > states;
  std::vector<size_t> effectiveK;

  size_t stateSize() const { return this->states.size(); }

  void enumerateOp();

  void makePairList();

  void makePairToEmission();

  void findUniqueState();

  void findEffectiveK();

  std::vector<int> makeTmpRow();

  std::vector<size_t> findWhichIsOne(std::vector<int> tmpOp);

  bool twoVectorsAreSame(std::vector<int> vec1, std::vector<int> vec2);

  std::vector<std::string> getIBDconfigureHeader();

};


class Hprior {
#ifdef UNITTEST
  friend class TestHprior;
  friend class TestMcmcMachinery;
  friend class TestIBDpath;
#endif

  friend class IBDpath;

  friend class McmcMachinery;

  friend class DEploidIO;

  Hprior() = default;
  ~Hprior() = default;

  void buildHprior(size_t kStrain, std::vector<double> &plaf);

  IBDconfiguration ibdConfig;
  size_t kStrain_;

  void setKstrain(const size_t setTo) { this->kStrain_ = setTo; }

  size_t kStrain() const { return this->kStrain_; }

  size_t nLoci_;

  void setnLoci(const size_t setTo) { this->nLoci_ = setTo; }

  size_t nLoci() const { return this->nLoci_; }

  std::vector<double> plaf_;
  std::vector<std::vector<double> > priorProb; // size: nState x nLoci
  std::vector<std::vector<double> > priorProbTrans; // size: nLoci x nState
  void transposePriorProbs();

  std::vector<size_t> stateIdx; // size: nState
  std::vector<size_t> stateIdxFreq;

  std::vector<std::vector<int> > hSet; // size: nState x kStrain
  size_t nState_;

  size_t nState() const { return this->nState_; }

  std::vector<size_t> effectiveK;

  size_t nPattern() const { return this->effectiveK.size(); }

  std::vector<std::string> getIBDconfigureHeader();
};


class IBDpath {
#ifdef UNITTEST
  friend class TestIBDpath;
#endif

  friend class McmcMachinery;

  friend class DEploidIO;

  std::shared_ptr<RandomGenerator> ibdRg_;

  double fSum;
  Hprior hprior;
  IBDrecombProbs ibdRecombProbs;
  std::vector<std::vector<double> > ibdTransProbs;
  std::vector<std::vector<double> > fm;
  std::vector<double> fSumState;
  std::vector<size_t> ibdConfigurePath;

  std::vector<std::vector<double> > bwd;
  std::vector<std::vector<double> > fwdbwd;

  IBDpath() = default;
  ~IBDpath() = default;

  size_t kStrain_;

  void setKstrain(const size_t setTo) { this->kStrain_ = setTo; }

  size_t kStrain() const { return this->kStrain_; }

  size_t nLoci_;

  void setNLoci(const size_t setTo) { this->nLoci_ = setTo; }

  size_t nLoci() const { return this->nLoci_; }

  double theta_;

  void setTheta(const double setTo) { this->theta_ = setTo; }

  double theta() const { return this->theta_; }

  std::vector<double> currentIBDpathChangeAt;

  std::vector<std::vector<double> > llkSurf;
  std::vector<int> uniqueEffectiveKCount;

  std::vector<double> IBDpathChangeAt;

  // Methods
  void computeAndUpdateTheta();

  void updateFmAtSiteI(std::vector<double> &prior,
                       std::vector<double> &llk);

  void ibdSamplePath(std::vector<double> statePrior);

  void makeIbdTransProbs();

  std::vector<double> computeEffectiveKPrior(double theta);

  std::vector<double> computeStatePrior(std::vector<double> effectiveKPrior);

  void makeLlkSurf(std::vector<double> altCount,
                   std::vector<double> refCount,
                   double scalingConst = 100.0,
                   double err = 0.01,
                   size_t gridSize = 99);

  void computeUniqueEffectiveKCount();

  std::vector<double> computeLlkOfStatesAtSiteI(std::vector<double> proportion, size_t siteI, double err = 0.01);

  std::vector<size_t> findWhichIsSomething(std::vector<size_t> tmpOp, size_t something);

  // For painting IBD
  void buildPathProbabilityForPainting(std::vector<double> proportion);

  void computeIbdPathFwdProb(std::vector<double> proportion, std::vector<double> statePrior);

  void computeIbdPathBwdProb(std::vector<double> proportion, std::vector<double> effectiveKPrior, std::vector<double> statePrior);

  void combineFwdBwd(std::vector<std::vector<double>> &reshapedFwd, std::vector<std::vector<double>> &reshapedBwd);

  std::vector<std::vector<double> > reshapeProbs(std::vector<std::vector<double> > &probs);

  double bestPath(std::vector<double> proportion, double err = 0.01);

public:
  std::vector<std::string> getIBDprobsHeader();

  void init(DEploidIO &dEploidIO, std::shared_ptr<RandomGenerator> randomGenerator);
};




}   // organization level namespace
}   // project level namespace



#endif
