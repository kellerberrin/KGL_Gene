//
// Created by kellerberrin on 20/7/21.
//

#include "kgl_variant_sort_analysis.h"


namespace kgl = kellerberrin::genome;



kgl::EnsemblIndexMap  kgl::SortedVariantAnalysis::filterEnsembl(const std::vector<std::string>& ensembl_list) const {

  EnsemblIndexMap filtered_map;

  for (auto const& ensembl_code : ensembl_list) {

    auto lower_bound = ensembl_index_map_->lower_bound(ensembl_code);
    auto upper_bound = ensembl_index_map_->upper_bound(ensembl_code);

    while(lower_bound != upper_bound) {

      filtered_map.insert(*lower_bound);

      ++lower_bound;

    }


  }

  return filtered_map;

}
