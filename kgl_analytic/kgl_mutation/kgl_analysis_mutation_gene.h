//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_H
#define KGL_ANALYSIS_MUTATION_GENE_H

#include "kgl_genome_genome.h"
#include "kgl_ped_parser.h"
#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_stats.h"
#include "kgl_analysis_mutation_gene_clinvar.h"
#include "kgl_analysis_mutation_gene_variant.h"


namespace kellerberrin::genome {   //  organization::project level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneMutation {

public:

  GeneMutation() = default;
  ~GeneMutation() = default;

  GeneMutation(const GeneMutation&) = default;
  GeneMutation& operator=(const GeneMutation&) = default;



  GeneCharacteristic gene_characteristic;
  GeneVariants gene_variants;
  GeneClinvar clinvar;

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using EnsemblIndexMap = std::multimap<std::string, std::shared_ptr<const Variant>>;
// How variants are allocated to genes.
enum class VariantGeneMembership { BY_SPAN, BY_EXON, BY_ENSEMBL };

class GenomeMutation {

public:

  GenomeMutation() { gene_membership_ = VariantGeneMembership::BY_EXON; }
  ~GenomeMutation() = default;

  // This analysis is performed first
  bool genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_reference,
                       const std::shared_ptr<const GenomePEDData>& ped_data);

  // Then this analysis.
  bool variantAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                        const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                        const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                        const std::shared_ptr<const GenomePEDData>& ped_data);
  // Output to file.
  bool writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data, const std::string& out_file, char output_delimiter) const;

  static EnsemblIndexMap ensemblIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr);

private:

  std::vector<GeneMutation> gene_vector_;
  VariantGeneMembership gene_membership_;


  constexpr static const char* CONCAT_TOKEN = "&";
  constexpr static const char* VEP_ENSEMBL_FIELD_ = "Gene";


  static void writeHeader(const std::shared_ptr<const GenomePEDData>& ped_data,
                          std::ostream& out_file,
                          char output_delimiter,
                          const GeneMutation& gene_mutation);

  static std::shared_ptr<const ContigDB> getGeneSpan( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                      const GeneCharacteristic& gene_char);


  static std::shared_ptr<const ContigDB> getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                     const GeneCharacteristic& gene_char);

  static std::shared_ptr<const ContigDB> getGeneEnsembl( const EnsemblIndexMap& ensembl_index_map,
                                                         const GeneCharacteristic& gene_char);

  GeneMutation geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                 const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                 const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                 const std::shared_ptr<const GenomePEDData>& ped_data,
                                 const EnsemblIndexMap& ensembl_index_map,
                                 GeneMutation gene_mutation);



};




} // namespace




#endif //KGL_ANALYSIS_MUTATION_GENE_H
