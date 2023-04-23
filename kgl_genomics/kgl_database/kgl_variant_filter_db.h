//
// Created by kellerberrin on 22/04/23.
//

#ifndef KGL_VARIANT_FILTER_DB_H
#define KGL_VARIANT_FILTER_DB_H


#include "kgl_variant_db_genome.h"
#include "kel_utility.h"


namespace kellerberrin::genome {   //  organization::project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter a population for genomes.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomeFilter : public FilterGenomes {

public:

  explicit GenomeFilter(const std::vector<GenomeId_t> &filtered_genomes) {

    filterName("GenomeFilter");

    for (auto const &genome_id: filtered_genomes) {

      filter_genomes_.insert(genome_id);

    }

  }
  ~GenomeFilter() override = default;

  [[nodiscard]] bool applyFilter(const GenomeDB& genome) const override { return filter_genomes_.contains(genome.genomeId()); }

  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<GenomeFilter>(*this); }

private:

  std::unordered_set<GenomeId_t> filter_genomes_;

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Filter for contigs within a genome.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFilter : public FilterContigs {

public:

  explicit ContigFilter(const std::vector<ContigId_t> &filtered_contigs) {

    filterName("Contig Filter");

    for (auto const &contig_id: filtered_contigs) {

      filter_contigs_.insert(contig_id);

    }

  }
  ~ContigFilter() override = default;

  [[nodiscard]] bool applyFilter(const ContigDB& contig) const override { return filter_contigs_.contains(contig.contigId()); }

  [[nodiscard]] std::shared_ptr<BaseFilter> clone() const override { return std::make_shared<ContigFilter>(*this); }

private:

  std::unordered_set<ContigId_t> filter_contigs_;

};






} // End namespace.

#endif //KGL_VARIANT_FILTER_DB_H
