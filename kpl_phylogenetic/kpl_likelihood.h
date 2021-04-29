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


namespace kellerberrin::phylogenetic {   //  organization level namespace

class Likelihood {

public:

  using SharedPtr = std::shared_ptr<Likelihood>;

  Likelihood();
  ~Likelihood();

  void setRooted(bool is_rooted);
  void setPreferGPU(bool prefer_gpu);
  void setAmbiguityEqualsMissing(bool ambig_equals_missing);

  bool usingStoredData() const { return _using_data; }
  void useStoredData(bool using_data) { _using_data = using_data; }

  void useUnderflowScaling(bool do_scaling);

  std::string beagleLibVersion() const;
  std::string availableResources() const;
  std::string usedResources() const;

  void initBeagleLib();
  void finalizeBeagleLib(bool use_exceptions);

  double calcLogLikelihood(const Tree& tree);

  Data::SharedPtr getData();
  void setData(Data::SharedPtr d);

  void clear();

  unsigned calcNumEdgesInFullyResolvedTree() const;
  unsigned calcNumInternalsInFullyResolvedTree() const;

  std::shared_ptr<Model> getModel();
  void setModel(std::shared_ptr<Model> model_ptr);

private:

  struct InstanceInfo {
    int handle{-1};
    int resourcenumber{-1};
    std::string resourcename;
    unsigned nstates{0};
    unsigned nratecateg{0};
    unsigned npatterns{0};
    unsigned partial_offset{0};
    unsigned tmatrix_offset{0};
    bool invarmodel{false};
    std::vector<unsigned> subsets;
  };

  using instance_pair_t = std::pair<unsigned, int>;

  unsigned getPartialIndex(Node::PtrNode  nd, InstanceInfo & info) const;
  unsigned getTMatrixIndex(Node::PtrNode  nd, InstanceInfo & info, unsigned subset_index) const;
  void updateInstanceMap(instance_pair_t & p, unsigned subset);
  void newInstance(unsigned nstates, int nrates, std::vector<unsigned> & subset_indices);
  void setTipStates();
  void setTipPartials();
  void setPatternPartitionAssignments();
  void setPatternWeights();
  void setAmongSiteRateHeterogenetity();
  void setModelRateMatrix();
  void addOperation(InstanceInfo & info, Node::PtrNode  nd, Node::PtrNode  lchild, Node::PtrNode  rchild, unsigned subset_index);
  void defineOperations(const Tree& tree);
  void updateTransitionMatrices();
  void calculatePartials();
  double calcInstanceLogLikelihood(InstanceInfo & inst, const Tree& tree);

// Setup Beagle.
  [[nodiscard]] static int setBeagleEigenDecomposition(std::shared_ptr<const Model> model_ptr, int beagle_instance, unsigned subset, unsigned instance_subset);
  [[nodiscard]] static int setBeagleStateFrequencies(std::shared_ptr<const Model> model_ptr, int beagle_instance, unsigned subset, unsigned instance_subset);
  [[nodiscard]] static int setBeagleAmongSiteRateVariationRates(std::shared_ptr<const Model> model_ptr, int beagle_instance, unsigned subset, unsigned instance_subset);
  [[nodiscard]] static int setBeagleAmongSiteRateVariationProbs(std::shared_ptr<const Model> model_ptr, int beagle_instance, unsigned subset, unsigned instance_subset);


  std::shared_ptr<Model> _model;
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

};



} // end Namespace


#endif // KPL_LIKELIHOOD_H
