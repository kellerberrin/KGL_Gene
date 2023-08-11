//
// Created by kellerberrin on 22/04/23.
//

#ifndef KGL_VARIANT_FILTER_DB_H
#define KGL_VARIANT_FILTER_DB_H

#include "kgl_variant_filter_type.h"
#include "kgl_variant_db_genome.h"
#include "kel_utility.h"



namespace kellerberrin::genome {   //  organization::project level namespace



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filters a population to the listed genomes (if they exist).
// Note that this is a shallow copy of the original population.
// Use selfFilter() or deepCopy() to create a permanent population view.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GenomeListFilter : public FilterPopulations {

public:

  explicit GenomeListFilter(const std::vector<GenomeId_t>& genome_list) : genome_set_(genome_list.begin(), genome_list.end()) {}
  explicit GenomeListFilter(std::set<GenomeId_t> genome_set) : genome_set_(std::move(genome_set)) {}
  ~GenomeListFilter() override = default;

  [[nodiscard]] std::unique_ptr<PopulationDB> applyFilter(const PopulationDB& population) const override;
  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<GenomeListFilter>(genome_set_); }

private:

  std::set<GenomeId_t> genome_set_;

};




} // End namespace.



#endif //KGL_VARIANT_FILTER_DB_H
