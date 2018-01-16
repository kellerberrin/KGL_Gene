//
// Created by kellerberrin on 16/01/18.
//

#ifndef KGL_PHYLOGENETIC_GENE_H
#define KGL_PHYLOGENETIC_GENE_H



#include <memory>
#include <fstream>
#include "kgl_library/kgl_patterns.h"
#include "kgl_variant_compound.h"
#include "kgl_library/kgl_variant_db.h"
#include "kgl_library/kgl_filter.h"
#include "kgl_library/kgl_gff_fasta.h"
#include "kgl_statistics.h"



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This object implements a series of high level functions for the kgl_phylogenetic application.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class GeneAnalysis {

public:

  GeneAnalysis() = default;
  virtual ~GeneAnalysis() = default;

  static bool mutateGene(const ContigId_t& contig,
                         const FeatureIdent_t& gene,
                         const FeatureIdent_t& sequence,
                         std::shared_ptr<const PopulationVariant> population_ptr,
                         std::shared_ptr<const GenomeDatabase> genome_db_ptr);

  static bool mutateGenomeGene(const GenomeId_t& genome,
                               const ContigId_t& contig,
                               const FeatureIdent_t& gene,
                               const FeatureIdent_t& sequence,
                               std::shared_ptr<const PopulationVariant> population_ptr,
                               std::shared_ptr<const GenomeDatabase> genome_db_ptr);


  static bool mutateRegion(const ContigId_t& contig,
                           ContigOffset_t offset,
                           ContigSize_t region_size,
                           std::shared_ptr<const PopulationVariant> population_ptr,
                           std::shared_ptr<const GenomeDatabase> genome_db_ptr);

  static bool mutateGenomeRegion(const GenomeId_t& genome,
                                 const ContigId_t& contig,
                                 const ContigOffset_t offset,
                                 const ContigSize_t region_size,
                                 std::shared_ptr<const PopulationVariant> population_ptr,
                                 std::shared_ptr<const GenomeDatabase> genome_db_ptr);

private:

  static bool mutateGenomeGene(const ContigId_t& contig,
                               const FeatureIdent_t& gene,
                               const FeatureIdent_t& sequence,
                               std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                               std::shared_ptr<const GenomeDatabase> genome_db_ptr);


  static bool mutateGenomeRegion(const ContigId_t& contig,
                                 const ContigOffset_t offset,
                                 const ContigSize_t region_size,
                                 std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                 std::shared_ptr<const GenomeDatabase> genome_db_ptr);

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_PHYLOGENETIC_GENE_H
