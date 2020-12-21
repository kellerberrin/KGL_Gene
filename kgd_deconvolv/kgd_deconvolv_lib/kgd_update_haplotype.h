

#ifndef KGD_HAP_H
#define KGD_HAP_H


#include <vector>
#include <iostream>
#include "kgd_utility.h"
#include "kgd_panel.h"


namespace kellerberrin {    // organization level namespace
namespace deconvolv {          // project level namespace



class UpdateHap {
#ifdef UNITTEST
  friend class TestUpdatePairHap;
  friend class TestUpdateSingleHap;
  friend class TestUpdateHap;
#endif


public:

  UpdateHap(size_t kStrain,
            size_t segmentStartIndex,
            size_t nLoci,
            std::shared_ptr<Panel> panel,
            double missCopyProb,
            double scalingFactor);

  virtual ~UpdateHap() = default;

  // Access functions.
  size_t nPanel() const { return nPanel_; }
  double getNewLLKIndex(size_t index) const { return newLLK[index]; }

  // Modification functions.
  void setPanelSize(const size_t setTo) { nPanel_ = setTo; }

protected:

  std::shared_ptr<Panel> panel_;
  double missCopyProb_;
  size_t kStrain_;
  size_t nPanel_;
  std::vector<double> newLLK;
  size_t segmentStartIndex_;
  size_t nLoci_;
  std::vector<std::vector<double> > emission_;
  double scalingFactor_;

  double scalingFactor() const { return scalingFactor_; }
  void setScalingFactor(const double setTo) { scalingFactor_ = setTo; }


};




}   // organization level namespace
}   // project level namespace



#endif
