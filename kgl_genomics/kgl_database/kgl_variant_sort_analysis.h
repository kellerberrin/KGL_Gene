//
// Created by kellerberrin on 20/7/21.
//

#ifndef KGL_VARIANT_SORT_ANALYSIS_H
#define KGL_VARIANT_SORT_ANALYSIS_H


#include "kgl_variant_sort.h"


namespace kellerberrin::genome {   //  organization::project level namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Sorts variants by Ensembl and "rsXXXXXX" Ids for further analysis.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class SortedVariantAnalysis {

public:

  explicit SortedVariantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr)
    : ensembl_index_map_(VariantSort::ensemblIndex(population_ptr)) {}
  ~SortedVariantAnalysis() = default;

  // Access the ensembl map.
  [[nodiscard]] const std::shared_ptr<const EnsemblIndexMap>& ensemblMap() const { return ensembl_index_map_; }
  // Filter on on list of Ensembl Codes.
  [[nodiscard]] EnsemblIndexMap filterEnsembl(const std::vector<std::string>& ensembl_list) const;

private:

  // A population of variants indexed by Ensembl gene code from the vep field.
  const std::shared_ptr<const EnsemblIndexMap> ensembl_index_map_;


};



} // namespace




#endif //KGL_VARIANT_SORT_ANALYSIS_H
