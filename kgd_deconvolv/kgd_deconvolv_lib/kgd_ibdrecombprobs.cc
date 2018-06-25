//
// Created by kellerberrin on 18/05/18.
//

#include <cmath>
#include "kgd_ibdrecombprobs.h"

namespace kgd = kellerberrin::deconvolv;


void kgd::IBDRecombProbs::computeRecombProbs(const std::vector<std::vector<int> >& position,
                                             size_t nLoci,
                                             double averageCentimorganDistance,
                                             double GFactor,
                                             bool useConstRecomb,
                                             double constRecombProb) {


  double averageMorganDistance = averageCentimorganDistance * 100;
  double geneticDistance;
  double rho;

  for (size_t i = 0; i < position.size(); i++) {

    for (size_t j = 1; j < position[i].size(); j++) {

      geneticDistance = (position[i][j] - position[i][j - 1]) / averageMorganDistance;
      rho = geneticDistance * GFactor;
      double pRecTmp = (useConstRecomb) ? constRecombProb : 1.0 - std::exp(-rho);
      recomb_probs_.push_back(pRecTmp);

    }

    recomb_probs_.push_back(1.0);

  }

  assert(recomb_probs_.size() == nLoci);

}


