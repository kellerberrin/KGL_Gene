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

#include "kgl_exec_env.h"
#include "kgd_update_haplotype.h"
#include <algorithm>    // std::reverse
#include <cstdlib>      // div


namespace kgd = kellerberrin::deploid;
namespace kgl = kellerberrin::genome;


kgd::UpdateHap::UpdateHap(const std::vector<double> &refCount,
                          const std::vector<double> &altCount,
                          const std::vector<double> &plaf,
                          const std::vector<double> &expectedWsaf,
                          const std::vector<double> &proportion,
                          const std::vector<std::vector<double> > &haplotypes,
                          std::shared_ptr<RandomGenerator> randomGenerator,
                          size_t segmentStartIndex,
                          size_t nLoci,
                          std::shared_ptr<Panel> panel,
                          double missCopyProb,
                          double scalingFactor) {

  panel_ = panel;

  if (panel_) {

    setPanelSize(panel_->truePanelSize());

  } else {

    setPanelSize(0);

  }

  kStrain_ = proportion.size();
  missCopyProb_ = missCopyProb;
  setScalingFactor(scalingFactor);
  recombRg_ = randomGenerator;
  recombLevel2Rg_ = randomGenerator;
  missCopyRg_ = randomGenerator;
  segmentStartIndex_ = segmentStartIndex;
  nLoci_ = nLoci;

}

void kgd::UpdateHap::core(std::vector<double> &refCount,
                     std::vector<double> &altCount,
                     std::vector<double> &plaf,
                     std::vector<double> &expectedWsaf,
                     std::vector<double> &proportion,
                     std::vector<std::vector<double> > &haplotypes) { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::calcExpectedWsaf(std::vector<double> &expectedWsaf,
                                      std::vector<double> &proportion,
                                      std::vector<std::vector<double> > &haplotypes) { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::calcHapLLKs(std::vector<double> &refCount,
                                 std::vector<double> &altCount) { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::buildEmission(double missCopyProb) { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::samplePaths() { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::addMissCopying(double missCopyProb) { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::updateLLK() { throw VirtualFunctionShouldNotBeCalled(); }

void kgd::UpdateHap::sampleHapIndependently(std::vector<double> &plaf) { throw VirtualFunctionShouldNotBeCalled(); }

