


#ifndef KGD_IBD_H
#define KGD_IBD_H


#include <vector>
#include <iostream>
#include <kgd_exceptions.h>
#include <sstream>
#include "kgd_utility.h"
#include "kgd_deploid_io.h"



namespace kellerberrin::deconvolv {          // project level namespace


// The IBDconfiguration is used for indexing.

class IBDconfiguration {

public:

  IBDconfiguration() = default;
  ~IBDconfiguration() = default;

  void buildIBDconfiguration(size_t k = 5);

  const std::vector<size_t>& effectiveK() const { return effectiveK_; }
  const std::vector<std::vector<size_t> >& states() const { return states_; }
  std::vector<std::string> getIBDconfigureHeader() const;

private:

  size_t kStrain_;
  std::vector<std::vector<size_t> > states_;
  std::vector<size_t> effectiveK_;


  void setKstrain(const size_t setTo) { kStrain_ = setTo; }

  size_t kStrain() const { return kStrain_; }

  void makePairList(std::vector< std::vector<size_t>>& pairs_list);

  void makePairToEmission(const std::vector< std::vector<size_t>>& pairs_list,
                          std::vector<std::vector<size_t> >& pairs_to_emission);

  void findEffectiveK();

  std::vector<size_t> makeEnumeratedArray();

  std::vector<size_t> activePairsArray(const std::vector<size_t>& pair_permute_row);

  static size_t nchoose2(size_t n);



};



}   // organization level namespace



#endif
