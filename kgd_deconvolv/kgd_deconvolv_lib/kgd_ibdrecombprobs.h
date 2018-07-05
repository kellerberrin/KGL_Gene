//
// Created by kellerberrin on 18/05/18.
//

#ifndef KGD_IBDRECOMBPROBS_H
#define KGD_IBDRECOMBPROBS_H

#include <vector>
#include "kgd_variant_index.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


/// Precompute (loci, loci-1) recombination probabilities for efficiency.
class IBDRecombProbs {

public:

  IBDRecombProbs(const std::vector<std::vector<size_t> >& position,
                 size_t nLoci,
                 double averageCentimorganDistance,
                 double Gfactor,
                 bool useConstRecomb,
                 double constRecombProb) {

    computeRecombProbs(position, nLoci, averageCentimorganDistance, Gfactor, useConstRecomb, constRecombProb);

  }

  ~IBDRecombProbs() = default;

  /// The recombination probability between loci and loci -1
  double getRecombProbAtLoci(size_t loci) const { return recomb_probs_[loci]; }
  /// The no recombination probability between loci and loci - 1
  double getNoRecombProbAtLoci(size_t loci) const { return (1.0 - recomb_probs_[loci]); }


private:

  std::vector<double> recomb_probs_;

  // Methods
  void computeRecombProbs(const std::vector<std::vector<size_t> >& position,
                          size_t nLoci,
                          double averageCentimorganDistance,
                          double Gfactor,
                          bool useConstRecomb,
                          double constRecombProb);

};



}   // organization level namespace
}   // project level namespace



#endif //KGD_IBDRECOMBPROBS_H
