//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_H
#define KGL_ANALYSIS_MUTATION_GENE_H

#include "kgl_genome_genome.h"
#include "kgl_Hsgenealogy_parser.h"
#include "kgl_uniprot_parser.h"
#include "kgl_variant_sort.h"
#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_stats.h"
#include "kgl_analysis_mutation_gene_clinvar.h"
#include "kgl_analysis_mutation_gene_variant.h"
#include "kgl_analysis_mutation_gene_ontology.h"


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
  OntologyStats ontology;

};




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// How variants are allocated to genes.
enum class VariantGeneMembership { BY_SPAN, BY_EXON, BY_ENSEMBL };
// By Span is all variants with in intron+exon span of the gene.
// By Exon is all variants within the exon region of the gene.
// By Ensembl looks up the Ensembl gene code in the variant vep.
// By EnsemblSummary is the same as the above but does not use the variant profile of a supplied population data file.
// Instead the statistics are generated directly from the summary (unphased single genome) data.

using EnsemblHashMap = std::unordered_map<std::string, const std::shared_ptr<const Variant>>;

class GenomeMutation {

public:

  explicit GenomeMutation(VariantGeneMembership gene_membership) : gene_membership_(gene_membership) {

    analysisType();


  }
  ~GenomeMutation() = default;

  // This analysis is performed first and only once.
  bool genomeAnalysis( const std::vector<std::string>& target_genes,
                       const std::shared_ptr<const GenomeReference>& genome_reference,
                       const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                       const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr,
                       const std::shared_ptr<const UniprotResource>& nomenclature_ptr);

  // Then this analysis.
  bool variantAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                        const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                        const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                        const std::shared_ptr<const HsGenomeAux>& genome_aux_data);

  // Finally, output to file.
  bool writeOutput(const std::shared_ptr<const HsGenomeAux>& genome_aux_data, const std::string& out_file, char output_delimiter) const;


private:

  std::vector<GeneMutation> gene_vector_;
  VariantGeneMembership gene_membership_;
  GeneEthnicitySex ethnic_statistics_;
  std::atomic<size_t> gene_variant_count_{0};
  std::atomic<size_t> ensembl_variant_count_{0};
  std::atomic<size_t> var_found_count_{0};
  std::atomic<size_t> var_checked_count_{0};


  static void writeHeader(const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                          std::ostream& out_file,
                          char output_delimiter,
                          const GeneMutation& gene_mutation);

  static std::shared_ptr<const ContigDB> getGeneSpan( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                      const GeneCharacteristic& gene_char);


  static std::shared_ptr<const ContigDB> getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                                     const GeneCharacteristic& gene_char);

  static std::shared_ptr<const ContigDB> getGeneEnsembl( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                         const EnsemblIndexMap& ensembl_index_map,
                                                         const GeneCharacteristic& gene_char);

  [[nodiscard]] std::shared_ptr<const ContigDB> getGeneEnsemblAlt( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                                   const EnsemblHashMap& ensembl_hash_map,
                                                                   const GeneCharacteristic& gene_char);

  GeneMutation geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                 const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                 const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                 const std::shared_ptr<const HsGenomeAux>& genome_aux_data,
                                 const std::shared_ptr<const EnsemblIndexMap>& ensembl_index_map_ptr,
                                 GeneMutation gene_mutation);

  void analysisType();

  // return order: hgnc_id, ensembl_id
  std::tuple<std::string, std::string> getNomenclature( const std::shared_ptr<const UniprotResource>& nomenclature_ptr,
                                                        const std::shared_ptr<const GeneFeature>& gene_ptr);

  // Set up Ensembl map.
  void getGeneEnsemblHashMap( const EnsemblIndexMap& ensembl_index_map,
                              const GeneCharacteristic& gene_char,
                              EnsemblHashMap& ensembl_hash_map,
                              ContigOffset_t& lower_bound,
                              ContigOffset_t& upper_bound);


  };




} // namespace




#endif //KGL_ANALYSIS_MUTATION_GENE_H
