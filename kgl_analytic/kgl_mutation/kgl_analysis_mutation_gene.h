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

// Variants indexed by their Ensembl gene id retrieved from the vep field.
using EnsemblIndexMap = std::multimap<std::string, std::shared_ptr<const Variant>>;
// How variants are allocated to genes.
enum class VariantGeneMembership { BY_SPAN, BY_EXON, BY_ENSEMBL };

class GenomeMutation {

public:

  explicit GenomeMutation(VariantGeneMembership gene_membership) : gene_membership_(gene_membership) { analysisType(); }
  ~GenomeMutation() = default;

  // This analysis is performed first and only once.
  bool genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_reference,
                       const std::shared_ptr<const GenomePEDData>& ped_data,
                       const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);

  // Then this analysis.
  bool variantAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                        const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                        const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                        const std::shared_ptr<const GenomePEDData>& ped_data,
                        const std::shared_ptr<const kol::OntologyDatabase>& ontology_db_ptr);

  // Finally, output to file.
  bool writeOutput(const std::shared_ptr<const GenomePEDData>& ped_data, const std::string& out_file, char output_delimiter) const;

  // A population of variants indexed by Ensembl gene code from the vep field. This may be useful elsewhere.
  [[nodiscard]] static std::shared_ptr<const EnsemblIndexMap> ensemblIndex(const std::shared_ptr<const PopulationDB>& unphased_population_ptr);

private:

  std::vector<GeneMutation> gene_vector_;
  bool analysis_initialized_{false};   // Only execute genomeAnalysis once
  VariantGeneMembership gene_membership_;

  constexpr static const char* CONCAT_TOKEN = "&";
  constexpr static const char* VEP_ENSEMBL_FIELD_ = "Gene";
  constexpr static const size_t PMR_BUFFER_SIZE_ = 4096;

  static void writeHeader(const std::shared_ptr<const GenomePEDData>& ped_data,
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

  static std::shared_ptr<const ContigDB> getGeneEnsemblSpan( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                             const EnsemblIndexMap& ensembl_index_map,
                                                             const GeneCharacteristic& gene_char);

  GeneMutation geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                 const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                 const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                 const std::shared_ptr<const GenomePEDData>& ped_data,
                                 const std::shared_ptr<const EnsemblIndexMap>& ensembl_index_map_ptr,
                                 GeneMutation gene_mutation);

  void analysisType();

};




} // namespace




#endif //KGL_ANALYSIS_MUTATION_GENE_H
