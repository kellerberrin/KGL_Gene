//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_PARTITION_H
#define KPL_PARTITION_H


#include "kpl_geneticcode.h"
#include "kpl_genetictype.h"
#include "kpl_xstrom.h"

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <tuple>
#include <limits>
#include <cmath>
#include <regex>



namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace


class Partition {

private:

  using regex_match_t = std::match_results<std::string::const_iterator>::const_reference;
  using subset_range_start_t = unsigned;
  using subset_range_end_t = unsigned;
  using subset_range_step_t = unsigned;
  using subset_index_t = unsigned;
  using subset_range_t = std::tuple<subset_range_start_t, subset_range_end_t, subset_range_step_t, subset_index_t>;
  using datatype_vect_t = std::vector<DataType> ;
  using subset_sizes_vect_t = std::vector<unsigned>;
  using subset_names_vect_t = std::vector<std::string>;

public:

  using SharedPtr = std::shared_ptr<Partition>;
  using partition_t = std::vector<subset_range_t>;

  Partition();
  ~Partition() = default;

  // Access
  [[nodiscard]] unsigned getNumSites() const;
  [[nodiscard]] unsigned getNumSubsets() const;
  [[nodiscard]] std::string getSubsetName(unsigned subset) const;

  [[nodiscard]] const partition_t& getSubsetRangeVect() const;

  [[nodiscard]] unsigned findSubsetByName(const std::string & subset_name) const;
  [[nodiscard]] unsigned findSubsetForSite(unsigned site_index) const;
  [[nodiscard]] bool siteInSubset(unsigned site_index, unsigned subset_index) const;
  [[nodiscard]] DataType getDataTypeForSubset(unsigned subset_index) const;
  [[nodiscard]] const datatype_vect_t& getSubsetDataTypes() const;

  [[nodiscard]] unsigned numSitesInSubset(unsigned subset_index) const;
  [[nodiscard]] subset_sizes_vect_t calcSubsetSizes() const;

  // Modify
  void defaultPartition(unsigned nsites = _INFINITY);
  void parseSubsetDefinition(std::string & s);
  void finalize(unsigned nsites);

private:

  unsigned _num_sites;
  unsigned _num_subsets;
  subset_names_vect_t _subset_names;
  partition_t _subset_ranges;
  datatype_vect_t _subset_data_types;

  static constexpr const unsigned _INFINITY = std::numeric_limits<unsigned>::max();

  void clear();
  int  extractIntFromRegexMatch(regex_match_t s, unsigned min_value);
  void addSubsetRange(unsigned subset_index, std::string range_definition);
  void addSubset(unsigned subset_index, std::string subset_definition);

};



} // phylogenetic
} // kellerberrin


#endif // KPL_PARTITION_H
