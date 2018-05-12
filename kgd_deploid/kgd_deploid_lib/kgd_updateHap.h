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




}   // organization level namespace
}   // project level namespace



#endif
