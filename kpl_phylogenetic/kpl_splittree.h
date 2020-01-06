//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_SPLITTREE_H
#define KPL_SPLITTREE_H

#include <vector>
#include <memory>
#include <set>
#include <map>
#include <climits>
#include <cassert>


namespace kellerberrin { //  organization level namespace
namespace phylogenetic { // project level namespace


class Split {

public:

  Split();

  Split(const Split &other);

  ~Split();

  Split &operator=(const Split &other);

  bool operator==(const Split &other) const;

  bool operator<(const Split &other) const;

  void clear();

  void resize(unsigned nleaves);

  typedef unsigned long split_unit_t;
  typedef std::vector<split_unit_t> split_t;
  typedef std::set<Split> treeid_t;
  typedef std::map<treeid_t, std::vector<unsigned> > treemap_t;
  typedef std::tuple<unsigned, unsigned, unsigned> split_metrics_t;

  split_unit_t getBits(unsigned unit_index) const;

  bool getBitAt(unsigned leaf_index) const;

  void setBitAt(unsigned leaf_index);

  void addSplit(const Split &other);

  bool isEquivalent(const Split &other) const;

  bool isCompatible(const Split &other) const;

  bool conflictsWith(const Split &other) const;

  std::string createPatternRepresentation() const;

  split_metrics_t getSplitMetrics() const;

private:

  split_unit_t _mask;
  split_t _bits;
  unsigned _bits_per_unit;
  unsigned _nleaves;

public:

  typedef std::shared_ptr<Split> SharedPtr;

};


} // phylogenetic
} // kellerberrin


#endif // KPL_SPLITTREE_H
