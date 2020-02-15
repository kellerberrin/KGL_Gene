//
// Created by kellerberrin on 18/12/19.
//

#ifndef KPL_MCMC_TREEUPDATER_H
#define KPL_MCMC_TREEUPDATER_H


#include "kpl_mcmc_updater.h"



namespace kellerberrin::phylogenetic {   //  organization level namespace


class Chain;

class TreeUpdater : public Updater {

//  friend class Chain;

public:

  typedef std::shared_ptr< TreeUpdater > SharedPtr;

  TreeUpdater();
  ~TreeUpdater();

  double calcLogPrior() override;

private:

  virtual void revert();
  virtual void proposeNewState();

  Node::PtrNode chooseRandomChild(Node::PtrNode  x, Node::PtrNode  avoid, bool parent_included);

  virtual void reset();

  void starTreeMove();

  double _orig_edgelen_top;
  double _orig_edgelen_middle;
  double _orig_edgelen_bottom;

  unsigned _case;
  bool _topology_changed;
  Node::PtrNode _x;
  Node::PtrNode _y;
  Node::PtrNode _a;
  Node::PtrNode _b;

  bool _star_tree_move;

};


} // end namespace


#endif //KPL_MCMC_TREEUPDATER_H
