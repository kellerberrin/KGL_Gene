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

#ifndef KGD_MCMC_H
#define KGD_MCMC_H


#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "kgd_random_generator.h"
#include "kgd_deploid_io.h"
#include "kgd_panel.h"
#include "randomSample.hpp"   // src/codeCogs/randomSample.hpp
#include "kgd_ibdpath.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



class McmcSample {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif

public:

  McmcSample() = default;
  ~McmcSample() = default;

  void clear() {

    proportion.clear();
    sumLLKs.clear();
    moves.clear();

  }

  std::vector<double> siteOfTwoSwitchOne;
  std::vector<double> siteOfTwoMissCopyOne;
  std::vector<double> siteOfTwoSwitchTwo;
  std::vector<double> siteOfTwoMissCopyTwo;
  std::vector<double> siteOfOneSwitchOne;
  std::vector<double> siteOfOneMissCopyOne;

  std::vector<double> currentsiteOfTwoSwitchOne;
  std::vector<double> currentsiteOfTwoMissCopyOne;
  std::vector<double> currentsiteOfTwoSwitchTwo;
  std::vector<double> currentsiteOfTwoMissCopyTwo;
  std::vector<double> currentsiteOfOneSwitchOne;
  std::vector<double> currentsiteOfOneMissCopyOne;
  std::vector<std::vector<double> > proportion;
  std::vector<std::vector<double> > hap;
  std::vector<double> sumLLKs;

  std::vector<int> moves;

};


class McmcMachinery {
#ifdef UNITTEST
  friend class TestMcmcMachinery;
#endif

public:

  McmcMachinery(std::shared_ptr<DEploidIO> dEplioidIO,
                std::shared_ptr<McmcSample> mcmcSample,
                std::shared_ptr<RandomGenerator> randomGenerator,
                bool useIBD = false);
  ~McmcMachinery() =default;

  void runMcmcChain(bool showProgress = true, bool useIBD = false, bool notInR = true);

private:

  std::shared_ptr<McmcSample> mcmcSample_;
  /* Variables */
  std::shared_ptr<DEploidIO> dEploidIO_;
  std::shared_ptr<Panel> panel_;
  size_t kStrain_;


  size_t nLoci_;

  double burnIn_;
  size_t maxIteration_;
  size_t mcmcThresh_;
  size_t McmcMachineryRate_;
  int eventInt_;

  size_t strainIndex_;
  size_t strainIndex1_;
  size_t strainIndex2_;

  size_t seed_;
  std::shared_ptr<RandomGenerator> hapRg_;
  std::shared_ptr<RandomGenerator> mcmcEventRg_;
  std::shared_ptr<RandomGenerator> propRg_;
  std::shared_ptr<RandomGenerator> initialHapRg_;

  //std::normal_distribution<double>* initialTitre_normal_distribution_;// (MN_LOG_TITRE, SD_LOG_TITRE);
  //std::normal_distribution<double>* deltaX_normal_distribution_;// (0, 1/PROP_SCALE);
  std::shared_ptr<StandNormalRandomSample> stdNorm_;

  double MN_LOG_TITRE;
  double SD_LOG_TITRE;
  double PROP_SCALE;

  size_t currentMcmcIteration_;
  std::vector<double> currentTitre_;
  double currentLogPriorTitre_;
  std::vector<double> currentProp_;
  std::vector<double> currentLLks_;
  std::vector<std::vector<double> > currentHap_;
  std::vector<double> currentExpectedWsaf_;
  std::vector<double> cumExpectedWsaf_;
  bool recordingMcmcBool_;

  IBDpath ibdPath; /* IBD */

  int acceptUpdate;

  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }


  void setNLoci(const size_t setTo) { nLoci_ = setTo; }

  size_t nLoci() const { return nLoci_; }

  double initialTitreNormalVariable() { return stdNorm_->genReal() * SD_LOG_TITRE + MN_LOG_TITRE; }

  //double deltaXnormalVariable(){ return stdNorm_->genReal() * 1.0/PROP_SCALE + MN_LOG_TITRE; }
  double deltaXnormalVariable() { return stdNorm_->genReal() * SD_LOG_TITRE * 1.0 / PROP_SCALE + MN_LOG_TITRE; }

  /* Methods */
  void calcMaxIteration(size_t nSample, size_t McmcMachineryRate, double burnIn);

  /* Initialize */
  void initializeMcmcChain(bool useIBD);

  void initializeProp();

  void initializeTitre();

  void initializeHap();

  void initializellk();

  void initializeExpectedWsaf();

  std::vector<double> calcExpectedWsaf(std::vector<double> &proportion);

  std::vector<double> titre2prop(std::vector<double> &tmpTitre);

  double calcLogPriorTitre(std::vector<double> &tmpTitre);

  double rBernoulli(double p);

  void printArray(std::vector<double> array) {

    for (auto const &value : array) {

      std::cout << value << " ";

    }

    std::cout << std::endl;

  }

  void sampleMcmcEvent(bool useIBD = false);

  void recordMcmcMachinery();


  void writeLastFwdProb(bool useIBD);

  void updateReferencePanel(size_t inbreedingPanelSizeSetTo, size_t excludedStrain);

  void initializeUpdateReferencePanel(size_t inbreedingPanelSizeSetTo);

  void computeDiagnostics();

  std::vector<double> computeLlkAtAllSites(double err = 0.01);

  std::vector<double> averageProportion(std::vector<std::vector<double> > &proportion);

  void ibdInitializeEssentials();

  void makeLlkSurf(std::vector<double> altCount,
                   std::vector<double> refCount,
                   double scalingConst = 100.0,
                   double err = 0.01,
                   size_t gridSize = 99);

  void ibdSampleMcmcEventStep();

  void initializePropIBD();

  void ibdSamplePath(std::vector<double> statePrior);

  void ibdUpdateHaplotypesFromPrior();

  void ibdUpdateProportionGivenHap(std::vector<double> &llkAtAllSites);
  //vector <double> getIBDprobsIntegrated(vector < vector <double> > &prob);

  /* Moves */
  void updateProportion();

  std::vector<double> calcTmpTitre();

  double deltaLLKs(std::vector<double> &newLLKs);

  void updateSingleHap();

  void findUpdatingStrainSingle();

  void updatePairHaps();

  //vector <size_t> sampleNoReplace(MersenneTwister* rg, vector <double> & proportion, size_t nSample );
  void findUpdatingStrainPair();

  /* Debug */
  bool doutProp();

  bool doutLLK();


};



}   // organization level namespace
}   // project level namespace



#endif

