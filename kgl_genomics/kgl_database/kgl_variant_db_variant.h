//
// Created by kellerberrin on 24/05/23.
//

#ifndef KGL_VARIANT_DB_VARIANT_H
#define KGL_VARIANT_DB_VARIANT_H

#include "kgl_variant_db_population.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


using VariantGenomeMap = std::map<GenomeId_t, uint32_t>;
class VariantDBRecord {

public:

  VariantDBRecord(const std::shared_ptr<const Variant>& variant_ptr, const VariantGenomeMap& genome_map) : variant_ptr_(variant_ptr), genome_map_(genome_map) {}
  ~VariantDBRecord() = default;

  [[nodiscard]] VariantGenomeMap& getMap() { return genome_map_; }
  [[nodiscard]] const Variant& getVariant() const { return *variant_ptr_; }

  [[nodiscard]] bool addGenomeCount(GenomeId_t genome, uint32_t variant_count);

private:

  std::shared_ptr<const Variant> variant_ptr_;
  VariantGenomeMap genome_map_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////

using VariantVariantMap = std::map<std::string, VariantDBRecord>; // Variants are indexed by their HGVS signature.
class VariantDBVariant {

public:

  VariantDBVariant(const std::shared_ptr<const PopulationDB>& population_ptr) { createVariantDB(population_ptr); }
  ~VariantDBVariant() = default;

  [[nodiscard]] const VariantVariantMap& getMap() const { return variant_map_; }

private:

  VariantVariantMap variant_map_;

  void createVariantDB(const std::shared_ptr<const PopulationDB>& population_ptr);

};



} // end namespace

#endif //KGL_VARIANT_DB_VARIANT_H
