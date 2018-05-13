//
// Created by kellerberrin on 12/05/18.
//

#ifndef KGD_UPDATEPAIRHAP_H
#define KGD_UPDATEPAIRHAP_H


#include <vector>
#include <iostream>
#include "kgd_utility.h"
#include "kgd_panel.h"
#include "kgd_update_haplotype.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace



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



#endif //KGL_KGD_UPDATEPAIRHAP_H
