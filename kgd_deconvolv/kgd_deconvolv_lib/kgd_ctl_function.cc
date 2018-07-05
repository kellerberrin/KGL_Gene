//
// Created by kellerberrin on 30/06/18.
//

#include "kgd_ctl_function.h"

namespace kgd = kellerberrin::deconvolv;


kgd::MixtureControlObj::MixtureControlObj() {


  setDoExportRecombProb(false);
  setrandomSeedWasGiven(false);
  setCompressVcf(false);
  setInitialPropWasGiven(false);
  setInitialHapWasGiven(false);

  setPleaseCheckInitialP(true);
  setExcludeSites(false);
  setUsePanel(true);

  setKStrainWasManuallySet(false);
  setKStrainWasSetByHap(false);
  setKStrainWasSetByProp(false);

  setDoUpdateProp(true);
  setDoUpdatePair(true);
  setDoUpdateSingle(true);
  setDoExportPostProb(false);
  setDoLsPainting(false);
  setDoIbdPainting(false);
  setUseIBD(false);
  setDoExportSwitchMissCopy(true);
  setDoAllowInbreeding(false);

  setForbidCopyFromSame(false);

  setUseVcf(false);
  setDoExportVcf(false);
  setDoComputeLLK(false);

  setDirectData(false);

}