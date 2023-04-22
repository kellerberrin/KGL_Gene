//
// Created by kellerberrin on 27/6/21.
//

#ifndef KGL_VARIANT_SORT_H
#define KGL_VARIANT_SORT_H


#include "kgl_variant.h"
#include "kgl_variant_db_population.h"

#include <map>
#include <string>
#include <memory>


namespace kellerberrin::genome {   //  organization::project level namespace

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convenience object sorts variants according to Ensembl Gene Id in the VEP field or Variant id.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


// filter indexed by their Ensembl gene id retrieved from the vep field.
using EnsemblIndexMap = std::multimap<std::string, std::shared_ptr<const Variant>>;

// filter indexed by their variant id ('rsXXXXXXXXX').
using VariantIdIndexMap = std::map<std::string, std::shared_ptr<const Variant>>;

// filter indexed by genome and then by variant id ('rsXXXXXXXXX').
using VariantGenomeIndexMap = std::map<std::string,  std::shared_ptr<VariantIdIndexMap>>;

// Indexed by variant_ID with value as Ensembl codes.
using VariantEnsemblIndexMap = std::map<std::string,  std::set<std::string>>;


// Static only functions, object cannot be created.
class VariantSort {

public:

  VariantSort()= delete;
  ~VariantSort() = delete;

  // A population of variants indexed by Ensembl gene code from the vep field.
  [[nodiscard]] static std::shared_ptr<EnsemblIndexMap> ensemblIndex(const std::shared_ptr<const PopulationDB>& population_ptr);

  // A population of variants indexed by Ensembl gene code from the vep field.
  // Only add the variants of the genes specified in the gene list.
  // An empty list adds all variants.
  static void ensemblAddIndex( const std::shared_ptr<const PopulationDB>& population_ptr,
                               const std::vector<std::string>& ensembl_gene_list,
                               std::shared_ptr<EnsemblIndexMap>& indexMap);

  [[nodiscard]] static size_t nonEnsemblIdentifiers(const EnsemblIndexMap& index_map);

  // A population of variants indexed by variant id.
  [[nodiscard]] static std::shared_ptr<VariantIdIndexMap> variantIdIndex(const std::shared_ptr<const PopulationDB>& population_ptr);

  // filter indexed by genome and then by variant id ('rsXXXXXXXXX').
  [[nodiscard]] static std::shared_ptr<VariantGenomeIndexMap> variantGenomeIndex(const std::shared_ptr<const PopulationDB>& population_ptr);

  // Multithreaded version indexes by variant id ('rsXXXXXXXXX') using a thread for each genome.
  [[nodiscard]] static std::shared_ptr<VariantGenomeIndexMap> variantGenomeIndexMT(const std::shared_ptr<const PopulationDB>& population_ptr);

private:

  constexpr static const char* VEP_ENSEMBL_FIELD_ = "Gene";
  constexpr static const size_t PMR_BUFFER_SIZE_ = 4096;
  constexpr static const char* ENSEMBL_PREFIX_ = "ENSG";

};



} // namespace

#endif //KGL_VARIANT_SORT_H
