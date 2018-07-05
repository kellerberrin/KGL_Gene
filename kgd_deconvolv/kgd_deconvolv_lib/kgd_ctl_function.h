//
// Created by kellerberrin on 30/06/18.
//

#ifndef KGD_CTL_FUNCTION_H
#define KGD_CTL_FUNCTION_H


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace


class MixtureControlObj {

public:

  explicit MixtureControlObj();
  ~MixtureControlObj() = default;

  MixtureControlObj& operator=(const MixtureControlObj& copy) = default;

  // Set Control
  void setInitialPropWasGiven(bool setTo) { initialPropWasGiven_ = setTo; }
  void setDoUpdateProp(bool setTo) { doUpdateProp_ = setTo; }
  void setInitialHapWasGiven(bool setTo) { initialHapWasGiven_ = setTo; }
  void setUseVcf(bool useVcf) { useVcf_ = useVcf; }
  void setUsePanel(bool setTo) { usePanel_ = setTo; }
  void setDoExportPostProb(bool setTo) { doExportPostProb_ = setTo; }
  void setDoExportSwitchMissCopy(bool setTo) { doExportSwitchMissCopy_ = setTo; }
  void setDoAllowInbreeding(bool setTo) { doAllowInbreeding_ = setTo; }
  void setDoIbdPainting(bool setTo) { doIbdPainting_ = setTo; }
  void setExcludeSites(bool exclude) { excludeSites_ = exclude; }
  void setKStrainWasManuallySet(bool setTo) { kStrainWasManuallySet_ = setTo; }
  void setrandomSeedWasGiven(bool random) { randomSeedWasGiven_ = random; }
  void setDoExportVcf(bool exportVcf) { doExportVcf_ = exportVcf; }
  void setUseIBD(bool setTo) { useIBD_ = setTo; }
  void setDirectData(bool setTo) { directData_ = setTo; }

  // Get Control
  bool doExportPostProb() const { return doExportPostProb_; }
  bool doAllowInbreeding() const { return doAllowInbreeding_; }
  bool initialPropWasGiven() const { return initialPropWasGiven_; }
  bool initialHapWasGiven() const { return initialHapWasGiven_; }
  bool doUpdateProp() const { return doUpdateProp_; }
  bool doUpdateSingle() const { return doUpdateSingle_; }
  bool doUpdatePair() const { return doUpdatePair_; }
  bool forbidCopyFromSame() const { return forbidCopyFromSame_; }
  bool useConstRecomb() const { return useConstRecomb_; }
  bool useIBD() const { return useIBD_; }
  bool useVcf() const { return useVcf_; }
  bool compressVcf() const { return compressVcf_; }
  bool doExportVcf() const { return doExportVcf_; }
  bool randomSeedWasGiven() const { return randomSeedWasGiven_; }
  bool doExportRecombProb() const { return doExportRecombProb_; }
  bool doLsPainting() const { return doLsPainting_; }
  bool doIbdPainting() const { return doIbdPainting_; }
  bool doComputeLLK() const { return doComputeLLK_; }
  bool excludeSites() const { return excludeSites_; }
  bool usePanel() const { return usePanel_; }
  bool directData() const { return directData_; }

private:

  bool randomSeedWasGiven_;
  bool initialPropWasGiven_;
  bool pleaseCheckInitialP_;
  bool initialHapWasGiven_;
  bool kStrainWasManuallySet_;
  bool kStrainWasSetByHap_;
  bool kStrainWasSetByProp_;
  bool useConstRecomb_;
  bool forbidCopyFromSame_;
  bool doUpdateProp_;
  bool doUpdatePair_;
  bool doUpdateSingle_;
  bool doExportPostProb_;
  bool doExportSwitchMissCopy_;
  bool doAllowInbreeding_;
  bool doLsPainting_;
  bool doIbdPainting_;
  bool useIBD_;
  bool usePanel_;  // Panel related
  bool useVcf_;
  bool doExportVcf_;
  bool compressVcf_;
  bool doExportRecombProb_;
  bool doComputeLLK_;
  bool excludeSites_;
  bool directData_;   // Data is injected directly from calling program (kgl).


  bool doExportSwitchMissCopy() const { return doExportSwitchMissCopy_; }
  bool pleaseCheckInitialP() const { return pleaseCheckInitialP_; }
  bool kStrainWasSetByHap() const { return kStrainWasSetByHap_; }
  bool kStrainWasManuallySet() const { return kStrainWasManuallySet_; }
  bool kStrainWasSetByProp() const { return kStrainWasSetByProp_; }

  // Getters and Setters
  void setDoUpdateSingle(bool setTo) { doUpdateSingle_ = setTo; }
  void setDoUpdatePair(bool setTo) { doUpdatePair_ = setTo; }
  void setDoLsPainting(bool setTo) { doLsPainting_ = setTo; }
  void setPleaseCheckInitialP(bool setTo) { pleaseCheckInitialP_ = setTo; }
  void setCompressVcf(bool compress) { compressVcf_ = compress; }
  void setDoComputeLLK(bool setTo) { doComputeLLK_ = setTo; }
  void setKStrainWasSetByHap(bool setTo) { kStrainWasSetByHap_ = setTo; }
  void setKStrainWasSetByProp(bool setTo) { kStrainWasSetByProp_ = setTo; }
  void setForbidCopyFromSame(bool forbid) { forbidCopyFromSame_ = forbid; }
  void setDoExportRecombProb(bool exportRecombProb) { doExportRecombProb_ = exportRecombProb; }

};


}   // organization level namespace
}   // project level namespace


#endif //KGD_CTL_FUNCTION_H
