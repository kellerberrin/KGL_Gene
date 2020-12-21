
#include "kgd_deconvolv_app.h"
#include "kgd_update_haplotype.h"
#include <cstdlib>      // div


namespace kgd = kellerberrin::deconvolv;


kgd::UpdateHap::UpdateHap(size_t kStrain,
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

  kStrain_ = kStrain;
  missCopyProb_ = missCopyProb;
  setScalingFactor(scalingFactor);
  segmentStartIndex_ = segmentStartIndex;
  nLoci_ = nLoci;

}

