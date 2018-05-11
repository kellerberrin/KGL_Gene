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


#ifndef KGD_HAP_H
#define KGD_HAP_H


#include <vector>
#include <iostream>
#include "kgd_utility.h"
#include "kgd_panel.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



class UpdateHap {
#ifdef UNITTEST
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
  friend class TestUpdateHap;
#endif

  friend class McmcMachinery;

  friend class UpdateSingleHap;

  friend class UpdatePairHap;

  friend class DEploidIO;

public:

  size_t nPanel() const { return this->nPanel_; }

private:

  UpdateHap(std::vector<double> &refCount,
            std::vector<double> &altCount,
            std::vector<double> &plaf,
            std::vector<double> &expectedWsaf,
            std::vector<double> &proportion,
            std::vector<std::vector<double> > &haplotypes,
            std::shared_ptr<RandomGenerator> randomGenerator,
            size_t segmentStartIndex,
            size_t nLoci,
            std::shared_ptr<Panel> panel,
            double missCopyProb,
            double scalingFactor);

  virtual ~UpdateHap() = default;

  std::shared_ptr<Panel> panel_;
  double missCopyProb_;
  std::shared_ptr<RandomGenerator> recombRg_;
  std::shared_ptr<RandomGenerator> recombLevel2Rg_;
  std::shared_ptr<RandomGenerator> missCopyRg_;

  size_t kStrain_;
  size_t nPanel_;

  void setPanelSize(const size_t setTo) { this->nPanel_ = setTo; }

  std::vector<double> newLLK;

  size_t segmentStartIndex_;
  size_t nLoci_;

  std::vector<std::vector<double> > emission_;
  double scalingFactor_;

  double scalingFactor() const { return this->scalingFactor_; }

  void setScalingFactor(const double setTo) { this->scalingFactor_ = setTo; }

  // Methods
  virtual void core(std::vector<double> &refCount,
                    std::vector<double> &altCount,
                    std::vector<double> &plaf,
                    std::vector<double> &expectedWsaf,
                    std::vector<double> &proportion,
                    std::vector<std::vector<double> > &haplotypes);

  virtual void calcExpectedWsaf(std::vector<double> &expectedWsaf, std::vector<double> &proportion,
                                std::vector<std::vector<double> > &haplotypes);

  virtual void calcHapLLKs(std::vector<double> &refCount, std::vector<double> &altCount);

  virtual void buildEmission(double missCopyProb);

  // calcFwdProbs() differ for class UpdateSingleHap and UpdatePairHap
  //virtual void calcFwdProbs();
  virtual void samplePaths();

  virtual void addMissCopying(double missCopyProb);

  virtual void updateLLK();

  virtual void sampleHapIndependently(std::vector<double> &plaf);
};


class UpdateSingleHap : public UpdateHap {
#ifdef UNITTEST
  friend class TestUpdateSingleHap;
#endif

  friend class McmcMachinery;

  friend class DEploidIO;
  //public:
private:
  UpdateSingleHap(std::vector<double> &refCount,
                  std::vector<double> &altCount,
                  std::vector<double> &plaf,
                  std::vector<double> &expectedWsaf,
                  std::vector<double> &proportion,
                  std::vector<std::vector<double> > &haplotypes,
                  std::shared_ptr<RandomGenerator> randomGenerator,
                  size_t segmentStartIndex,
                  size_t nLoci,
                  std::shared_ptr<Panel> panel, double missCopyProb,
                  double scalingFactor,
                  size_t strainIndex);

  ~UpdateSingleHap() = default;

  std::vector<double> siteOfOneSwitchOne;
  std::vector<double> siteOfOneMissCopyOne;
  std::vector<std::vector<double> > fwdProbs_;
  std::vector<std::vector<double> > bwdProbs_;
  std::vector<std::vector<double> > fwdBwdProbs_;

  size_t strainIndex_;
  std::vector<double> expectedWsaf0_;
  std::vector<double> expectedWsaf1_;
  std::vector<double> llk0_;
  std::vector<double> llk1_;

  std::vector<double> path_;
  std::vector<double> hap_;

  // Methods
  void core(std::vector<double> &refCount,
            std::vector<double> &altCount,
            std::vector<double> &plaf,
            std::vector<double> &expectedWsaf,
            std::vector<double> &proportion,
            std::vector<std::vector<double> > &haplotypes);

  void painting(std::vector<double> &refCount,
                std::vector<double> &altCount,
                std::vector<double> &expectedWsaf,
                std::vector<double> &proportion,
                std::vector<std::vector<double> > &haplotypes);

  void calcExpectedWsaf(std::vector<double> &expectedWsaf, std::vector<double> &proportion,
                        std::vector<std::vector<double> > &haplotypes);

  void calcHapLLKs(std::vector<double> &refCount, std::vector<double> &altCount);

  void buildEmission(double missCopyProb);

  void buildEmissionBasicVersion(double missCopyProb);

  void calcFwdProbs();

  void calcBwdProbs();

  void calcFwdBwdProbs();

  void samplePaths();

  void addMissCopying(double missCopyProb);

  void sampleHapIndependently(std::vector<double> &plaf);

  void updateLLK();
};


class UpdatePairHap : public UpdateHap {
#ifdef UNITTEST
  friend class TestUpdatePairHap;
#endif

  friend class McmcMachinery;

  friend class DEploidIO;

public:

  UpdatePairHap(std::vector<double> &refCount,
                std::vector<double> &altCount,
                std::vector<double> &plaf,
                std::vector<double> &expectedWsaf,
                std::vector<double> &proportion,
                std::vector<std::vector<double> > &haplotypes,
                std::shared_ptr<RandomGenerator> randomGenerator,
                size_t segmentStartIndex,
                size_t nLoci,
                std::shared_ptr<Panel> panel, double missCopyProb,
                double scalingFactor, bool forbidCopyFromSame,
                size_t strainIndex1,
                size_t strainIndex2);

  ~UpdatePairHap() = default;

private:
  std::vector<double> siteOfTwoSwitchOne;
  std::vector<double> siteOfTwoMissCopyOne;
  std::vector<double> siteOfTwoSwitchTwo;
  std::vector<double> siteOfTwoMissCopyTwo;
  std::vector<std::vector<std::vector<double> > > fwdProbs_;

  size_t strainIndex1_;
  size_t strainIndex2_;
  bool forbidCopyFromSame_;

  std::vector<double> expectedWsaf00_;
  std::vector<double> expectedWsaf01_;
  std::vector<double> expectedWsaf10_;
  std::vector<double> expectedWsaf11_;
  std::vector<double> llk00_;
  std::vector<double> llk01_;
  std::vector<double> llk10_;
  std::vector<double> llk11_;
  std::vector<double> path1_;
  std::vector<double> path2_;
  std::vector<double> hap1_;
  std::vector<double> hap2_;

  // Methods
  void core(std::vector<double> &refCount,
            std::vector<double> &altCount,
            std::vector<double> &plaf,
            std::vector<double> &expectedWsaf,
            std::vector<double> &proportion,
            std::vector<std::vector<double> > &haplotypes);

  void calcExpectedWsaf(std::vector<double> &expectedWsaf, std::vector<double> &proportion,
                        std::vector<std::vector<double> > &haplotypes);

  void calcHapLLKs(std::vector<double> &refCount, std::vector<double> &altCount);

  void buildEmission(double missCopyProb);

  void calcFwdProbs(bool forbidCopyFromSame);

  void samplePaths();

  void addMissCopying(double missCopyProb);

  void sampleHapIndependently(std::vector<double> &plaf);

  void updateLLK();

  // Own methods
  std::vector<double> computeRowMarginalDist(std::vector<std::vector<double> > &probDist);

  std::vector<double> computeColMarginalDist(std::vector<std::vector<double> > &probDist);

  std::vector<size_t> sampleMatrixIndex(std::vector<std::vector<double> > &probDist);

};



}   // organization level namespace
}   // project level namespace



#endif
