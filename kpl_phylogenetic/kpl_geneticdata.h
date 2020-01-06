//
// Created by kellerberrin on 11/12/19.
//

#ifndef KPL_GENETICDATA_H
#define KPL_GENETICDATA_H

#include "kpl_geneticcode.h"
#include "kpl_genetictype.h"
#include "kpl_partition.h"
#include "kpl_xstrom.h"

#include "ncl/nxsmultiformat.h"

#include <boost/algorithm/string/join.hpp>
#include <boost/format.hpp>

#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <numeric>
#include <limits>
#include <map>


namespace kellerberrin {   //  organization level namespace
namespace phylogenetic {   // project level namespace



class Data {
public:
  typedef std::vector<std::string> taxon_names_t;
  typedef unsigned long long state_t;
  typedef std::vector<state_t> pattern_vect_t;
  typedef std::vector<state_t> monomorphic_vect_t;
  typedef std::vector<int> partition_key_t;
  typedef std::map<pattern_vect_t, unsigned> pattern_map_t;
  typedef std::vector<pattern_vect_t> data_matrix_t;
  typedef std::vector<pattern_map_t> pattern_map_vect_t;
  typedef std::vector<double> pattern_counts_t;
  typedef std::vector<unsigned> subset_end_t;
  typedef std::vector<unsigned> npatterns_vect_t;
  typedef std::pair<unsigned, unsigned> begin_end_pair_t;
  typedef std::shared_ptr<Data> SharedPtr;

  Data();

  ~Data();

  Partition::SharedPtr getPartition();

  void setPartition(Partition::SharedPtr partition);

  void getDataFromFile(const std::string filename);

  unsigned getNumSubsets() const;

  std::string getSubsetName(unsigned subset) const;

  unsigned getNumTaxa() const;

  const taxon_names_t &getTaxonNames() const;

  unsigned getNumPatterns() const;

  npatterns_vect_t calcNumPatternsVect() const;

  unsigned getNumPatternsInSubset(unsigned subset) const;

  unsigned getNumStatesForSubset(unsigned subset) const;

  unsigned calcSeqLen() const;

  unsigned calcSeqLenInSubset(unsigned subset) const;

  const data_matrix_t &getDataMatrix() const;

  begin_end_pair_t getSubsetBeginEnd(unsigned subset) const;

  const pattern_counts_t &getPatternCounts() const;

  const monomorphic_vect_t &getMonomorphic() const;

  const partition_key_t &getPartitionKey() const;

  std::string createTaxaBlock() const;

  std::string createTranslateStatement() const;


  void clear();

private:

  unsigned storeTaxonNames(NxsTaxaBlock *taxaBlock, unsigned taxa_block_index);

  unsigned
  storeData(unsigned ntax, unsigned nchar, NxsCharactersBlock *charBlock, NxsCharactersBlock::DataTypesEnum datatype);

  unsigned buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets);

  void updatePatternMap(Data::pattern_vect_t &pattern, unsigned subset);

  void compressPatterns();

  Partition::SharedPtr _partition;
  pattern_counts_t _pattern_counts;
  monomorphic_vect_t _monomorphic;
  partition_key_t _partition_key;
  pattern_map_vect_t _pattern_map_vect;
  taxon_names_t _taxon_names;
  data_matrix_t _data_matrix;
  subset_end_t _subset_end;
};



} // phylogenetic
} // kellerberrin


#endif //KPL_GENETICDATA_H
