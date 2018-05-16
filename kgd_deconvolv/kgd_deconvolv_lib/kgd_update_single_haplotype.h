//
// Created by kellerberrin on 12/05/18.
//

#ifndef KGD_UPDATESINGLEHAP_H
#define KGD_UPDATESINGLEHAP_H


#include "kgd_update_haplotype.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


class UpdateSingleHap : public UpdateHap {
#ifdef UNITTEST
  friend class TestUpdateSingleHap;
#endif

public:

  UpdateSingleHap(std::shared_ptr<RandomGenerator> randomGenerator,
                  size_t segmentStartIndex,
                  size_t nLoci,
                  size_t kStrain,
                  std::shared_ptr<Panel> panel,
                  double missCopyProb,
                  double scalingFactor,
                  size_t strainIndex);

  ~UpdateSingleHap() = default;

  const std::vector<std::vector<double> >& getFwdProbs() const { return fwdProbs_; }
  const std::vector<std::vector<double> >& getFwdBwdProbs() const { return fwdBwdProbs_; }
  double getHapIndex(size_t index) const { return hap_[index]; }
  double getOneSwitchOneIndex(size_t index) const { return siteOfOneSwitchOne[index]; }
  double getOneMissCopyOneIndex(size_t index) const { return siteOfOneMissCopyOne[index]; }

  // Methods
  void core(const std::vector<double> &refCount,
            const std::vector<double> &altCount,
            const std::vector<double> &plaf,
            const std::vector<double> &expectedWsaf,
            const std::vector<double> &proportion,
            const std::vector<std::vector<double> > &haplotypes);

  void painting(std::vector<double> &refCount,
                std::vector<double> &altCount,
                std::vector<double> &expectedWsaf,
                std::vector<double> &proportion,
                std::vector<std::vector<double> > &haplotypes);

private:

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


  void calcExpectedWsaf(const std::vector<double> &expectedWsaf,
                        const std::vector<double> &proportion,
                        const std::vector<std::vector<double> > &haplotypes);

  void calcHapLLKs(const std::vector<double> &refCount, const std::vector<double> &altCount);

  void buildEmission(double missCopyProb);

  void buildEmissionBasicVersion(double missCopyProb);

  void calcFwdProbs();

  void calcBwdProbs();

  void calcFwdBwdProbs();

  void samplePaths();

  void addMissCopying(double missCopyProb);

  void sampleHapIndependently(const std::vector<double> &plaf);

  void updateLLK();
};



}   // organization level namespace
}   // project level namespace



#endif //KGL_KGD_UPDATESINGLEHAP_H
