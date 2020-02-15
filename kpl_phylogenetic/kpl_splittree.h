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


namespace kellerberrin::phylogenetic { //  organization level namespace


class Split {

public:

  using split_unit_t = unsigned long ;
  using split_t = std::vector<split_unit_t> ;
  using treeid_t = std::set<Split> ;
  using treemap_t = std::map<treeid_t, std::vector<unsigned> > ;
  using split_metrics_t = std::tuple<unsigned, unsigned, unsigned>;

  Split();

  Split(const Split &other) = default;
  ~Split() = default;

  Split& operator=(const Split &copy) = default;

  bool operator==(const Split &other) const;

  bool operator<(const Split &other) const;

  void clear();

  void resize(unsigned nleaves);


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


} // end namespace


#endif // KPL_SPLITTREE_H
