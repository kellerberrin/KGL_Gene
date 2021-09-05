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



const std::shared_ptr<const kgl::VariantEnsemblIndexMap>& kgl::SortedVariantAnalysis::alleleEnsemblMap() const {

  if (variant_ensembl_index_map_) {

    return variant_ensembl_index_map_;

  }

  VariantEnsemblIndexMap variant_ensembl_map;

  for (auto const& [ensembl_code, variant_ptr] : *ensembl_index_map_) {

    if (not variant_ptr->identifier().empty() and not ensembl_code.empty()) {

      auto result = variant_ensembl_map.find(variant_ptr->identifier());
      if (result == variant_ensembl_map.end()) {

        variant_ensembl_map.emplace(variant_ptr->identifier(), std::set<std::string>{ensembl_code});

      } else {

        auto& [rs_key, ensembl_set] = *result;
        ensembl_set.insert(ensembl_code);

      }

    }

  }

  variant_ensembl_index_map_ = std::make_shared<const VariantEnsemblIndexMap>(std::move(variant_ensembl_map));

  return variant_ensembl_index_map_;

}
