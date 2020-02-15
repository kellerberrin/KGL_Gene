//
// Created by kellerberrin on 20/12/19.
//

#ifndef KPL_MCMC_POLYTOMYUPDATER_H
#define KPL_MCMC_POLYTOMYUPDATER_H

#include "kpl_mcmc_updater.h"

#include <vector>
#include <map>
#include <memory>


namespace kellerberrin::phylogenetic {   //  organization level namespace

class PolytomyUpdater : public Updater {

//  friend class Chain;

public:

  typedef std::vector<double> _partition_vect_t;
  typedef std::map<unsigned, _partition_vect_t > _partition_map_t;
  typedef std::vector<Node::PtrNode > _polytomy_vect_t;
  typedef std::shared_ptr< PolytomyUpdater > SharedPtr;

  PolytomyUpdater();
  ~PolytomyUpdater() override;

  double calcLogPrior() override;

private:

  void revert() override;
  void proposeNewState() override;
  void reset() override;

  void proposeAddEdgeMove(Node::PtrNode  nd);
  void proposeDeleteEdgeMove(Node::PtrNode  nd);

  _partition_vect_t& computePolytomyDistribution(unsigned nspokes);
  void refreshPolytomies();

  _partition_map_t _poly_prob;
  _polytomy_vect_t _polytomies;

  Node::PtrNode _orig_par;
  Node::PtrNode _orig_lchild;

  bool _add_edge_proposed;
  double _new_edge_proportion;
  double _orig_edge_proportion;
  double _tree_length;
  double _phi;
  unsigned _polytomy_size;
  unsigned _num_polytomies;

};



} // end namespace


#endif //KPL_MCMC_POLYTOMYUPDATER_H
