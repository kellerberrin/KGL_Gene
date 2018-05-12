//
// Created by kellerberrin on 12/05/18.
//

#ifndef KGD_UPDATESINGLEHAP_H
#define KGD_UPDATESINGLEHAP_H


#include "kgd_updateHap.h"


namespace kellerberrin {    // organization level namespace
namespace deploid {          // project level namespace


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



}   // organization level namespace
}   // project level namespace



#endif //KGL_KGD_UPDATESINGLEHAP_H
