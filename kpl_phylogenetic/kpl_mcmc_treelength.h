//
// Created by kellerberrin on 18/12/19.
//

#ifndef KPL_MCMC_TREELENGTH_H
#define KPL_MCMC_TREELENGTH_H

#include "kpl_mcmc_updater.h"


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class TreeLengthUpdater : public Updater {

public:

  typedef std::shared_ptr< TreeLengthUpdater > SharedPtr;

  TreeLengthUpdater();
  ~TreeLengthUpdater();

  virtual void clear();
  virtual void proposeNewState();
  virtual void revert();

  double calcLogPrior() override;

  void  pullFromModel();
  void  pushToModel();

private:

  double _prev_point;
  double _curr_point;

};


} // phylogenetic
} // kellerberrin



#endif //KPL_MCMC_TREELENGTH_H
