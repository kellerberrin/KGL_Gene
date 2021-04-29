//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_TREESUMMARY_H
#define KPL_TREESUMMARY_H


#include "kpl_splittree.h"
#include "kpl_treemanip.h"
#include "kpl_xstrom.h"

#include "ncl/nxsmultiformat.h"

#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>


namespace kellerberrin::phylogenetic {   //  organization level namespace


class TreeSummary {
public:
  TreeSummary() = default;
  ~TreeSummary() =default;

  void readTreefile(const std::string filename, unsigned skip);

  void showSummary() const;

  typename std::shared_ptr<Tree> getTree(unsigned index);

  std::string getNewick(unsigned index);

  void clear();

private:

  Split::treemap_t _treeIDs;
  std::vector<std::string> _newicks;

public:

  typedef std::shared_ptr<TreeSummary> SharedPtr;
};



} // end namespace


#endif // KPL_TREESUMMARY_H
