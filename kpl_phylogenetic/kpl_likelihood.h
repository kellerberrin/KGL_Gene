//
// Created by kellerberrin on 12/12/19.
//

#ifndef KPL_LIKELIHOOD_H
#define KPL_LIKELIHOOD_H

#include "kpl_model.h"
#include "kpl_tree.h"
#include "kpl_geneticdata.h"
#include "kpl_xstrom.h"

#include "libhmsbeagle/beagle.h"

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <map>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Likelihood {
public:
  Likelihood();
  ~Likelihood();

  void                                    setRooted(bool is_rooted);
  void                                    setPreferGPU(bool prefer_gpu);
  void                                    setAmbiguityEqualsMissing(bool ambig_equals_missing);

  bool                                    usingStoredData() const;
  void                                    useStoredData(bool using_data);

  void                                    useUnderflowScaling(bool do_scaling);

  std::string                             beagleLibVersion() const;
  std::string                             availableResources() const;
  std::string                             usedResources() const;

  void                                    initBeagleLib();
  void                                    finalizeBeagleLib(bool use_exceptions);

  double                                  calcLogLikelihood(Tree::SharedPtr tree);

  Data::SharedPtr                         getData();
  void                                    setData(Data::SharedPtr d);

  void                                    clear();

  unsigned                                calcNumEdgesInFullyResolvedTree() const;
  unsigned                                calcNumInternalsInFullyResolvedTree() const;

  Model::SharedPtr                        getModel();
  void                                    setModel(Model::SharedPtr m);

private:

  struct InstanceInfo {
    int handle;
    int resourcenumber;
    std::string resourcename;
    unsigned nstates;
    unsigned nratecateg;
    unsigned npatterns;
    unsigned partial_offset;
    unsigned tmatrix_offset;
    bool invarmodel;
    std::vector<unsigned> subsets;

    InstanceInfo() : handle(-1), resourcenumber(-1), resourcename(""), nstates(0), nratecateg(0), npatterns(0), partial_offset(0), tmatrix_offset(0), invarmodel(false) {}

  };

  typedef std::pair<unsigned, int>        instance_pair_t;

  unsigned                                getPartialIndex(Node::PtrNode  nd, InstanceInfo & info) const;
  unsigned                                getTMatrixIndex(Node::PtrNode  nd, InstanceInfo & info, unsigned subset_index) const;
  void                                    updateInstanceMap(instance_pair_t & p, unsigned subset);
  void                                    newInstance(unsigned nstates, int nrates, std::vector<unsigned> & subset_indices);
  void                                    setTipStates();
  void                                    setTipPartials();
  void                                    setPatternPartitionAssignments();
  void                                    setPatternWeights();
  void                                    setAmongSiteRateHeterogenetity();
  void                                    setModelRateMatrix();
  void                                    addOperation(InstanceInfo & info, Node::PtrNode  nd, Node::PtrNode  lchild, Node::PtrNode  rchild, unsigned subset_index);
  void                                    defineOperations(Tree::SharedPtr tree);
  void                                    updateTransitionMatrices();
  void                                    calculatePartials();
  double                                  calcInstanceLogLikelihood(InstanceInfo & inst, Tree::SharedPtr tree);

  Model::SharedPtr _model;
  std::vector<InstanceInfo> _instances;
  std::map<int, std::string> _beagle_error;
  std::map<int, std::vector<int> > _operations;
  std::map<int, std::vector<int> >  _pmatrix_index;
  std::map<int, std::vector<double> > _edge_lengths;

  std::vector<int> _subset_indices;
  std::vector<int> _parent_indices;
  std::vector<int> _child_indices;
  std::vector<int> _tmatrix_indices;
  std::vector<int> _weights_indices;
  std::vector<int> _freqs_indices;
  std::vector<int> _scaling_indices;

  Data::SharedPtr _data;
  unsigned _ntaxa;
  bool _rooted;
  bool _prefer_gpu;
  bool _ambiguity_equals_missing;

  bool _underflow_scaling;

  bool _using_data;

public:
  typedef std::shared_ptr< Likelihood >   SharedPtr;
};



} // Namespace phylogenetic
} // kellerberrin


#endif // KPL_LIKELIHOOD_H
