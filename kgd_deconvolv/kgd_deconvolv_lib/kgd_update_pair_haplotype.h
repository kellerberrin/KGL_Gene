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
namespace deconvolv {          // project level namespace



class UpdatePairHap : public UpdateHap {
#ifdef UNITTEST
  friend class TestUpdatePairHap;
#endif


public:

  UpdatePairHap(size_t segmentStartIndex,
                size_t nLoci,
                size_t kStrain,
                std::shared_ptr<Panel> panel,
                double missCopyProb,
                double scalingFactor,
                bool forbidCopyFromSame,
                size_t strainIndex1,
                size_t strainIndex2);

  ~UpdatePairHap() = default;

  // Methods
  void core(const std::vector<double> &refCount,
            const std::vector<double> &altCount,
            const std::vector<double> &plaf,
            const std::vector<double> &expectedWsaf,
            const std::vector<double> &proportion,
            const std::vector<std::vector<double> > &haplotypes);

  double getHap1Index(size_t index) const { return hap1_[index]; }
  double getHap2Index(size_t index) const { return hap2_[index]; }
  double getTwoSwitchOneIndex(size_t index) const { return siteOfTwoSwitchOne_[index]; }
  double getTwoMissCopyOneIndex(size_t index) const { return siteOfTwoMissCopyOne_[index]; }
  double getTwoSwitchTwoIndex(size_t index) const { return siteOfTwoSwitchTwo_[index]; }
  double getTwoMissCopyTwoIndex(size_t index) const { return siteOfTwoMissCopyTwo_[index]; }

private:

  std::vector<double> siteOfTwoSwitchOne_;
  std::vector<double> siteOfTwoMissCopyOne_;
  std::vector<double> siteOfTwoSwitchTwo_;
  std::vector<double> siteOfTwoMissCopyTwo_;
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


  void calcExpectedWsaf(const std::vector<double> &expectedWsaf,
                        const std::vector<double> &proportion,
                        const std::vector<std::vector<double> > &haplotypes);

  void calcHapLLKs(const std::vector<double> &refCount, const std::vector<double> &altCount);

  void buildEmission(double missCopyProb);

  void calcFwdProbs(bool forbidCopyFromSame);

  void samplePaths();

  void addMissCopying(double missCopyProb);

  void sampleHapIndependently(const std::vector<double> &plaf);

  void updateLLK();

  // Own methods
  std::vector<double> computeRowMarginalDist(std::vector<std::vector<double> > &probDist);

  std::vector<double> computeColMarginalDist(std::vector<std::vector<double> > &probDist);

  std::vector<size_t> sampleMatrixIndex(std::vector<std::vector<double> > &probDist);

};



}   // organization level namespace
}   // project level namespace



#endif //KGL_KGD_UPDATEPAIRHAP_H
