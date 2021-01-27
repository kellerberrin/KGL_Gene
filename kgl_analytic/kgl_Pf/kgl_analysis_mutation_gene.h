//
// Created by kellerberrin on 15/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_H
#define KGL_ANALYSIS_MUTATION_GENE_H

#include "kgl_genome_genome.h"
#include "kgl_ped_parser.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct GeneMutation {

  GenomeId_t genome;
  ContigId_t contig;
  FeatureIdent_t gene_id;
  std::string gene_name;
  std::string description;
  std::string biotype;
  bool valid_protein{false};
  std::string gaf_id;
  ContigOffset_t gene_begin{0};
  ContigOffset_t gene_end{0};
  ContigSize_t gene_span{0};
  ContigSize_t exons{0};
  ContigSize_t nucleotides{0};
  std::string strand;
  size_t sequences{0};
  std::string seq_name;
  size_t attribute_size{0};
  size_t unique_variants{0};
  size_t span_variant_count{0};
  size_t variant_count{0};
  size_t male_phase{0};  // Variants from the male phased (B) chromosome.
  size_t female_phase{0};  // Variants from the female phased (A) chromosome.
  size_t male_lof{0};          // Loss of gene function in the (B) chromosome.
  size_t female_lof{0};        // Los of function in the female (A) chromosome.
  size_t hom_lof{0};          // Loss of gene function in both chromosomes.
  size_t clinvar_count{0};
  std::string clinvar_desc;
  size_t male_high_effect{0};          // High Variant Impact in the (B) chromosome.
  size_t female_high_effect{0};        // High Variant Impact in the (A) chromosome.
  size_t hom_high_effect{0};          // High Impact in both in both chromosomes.
  size_t EAS{0};
  size_t EUR{0};
  size_t AFR{0};
  size_t AMR{0};
  size_t SAS{0};
  size_t genome_count{0};   // Total number of genomes.
  size_t male_value{0};    // Males that have values for this gene.
  size_t female_value{0};  // Females that have values for this gene.
  size_t genome_variant{0};  // Number of genomes that contain variants for this gene.
  size_t homozygous{0};
  size_t heterozygous{0};
  double indel{0.0};
  double transition{0.0};
  double transversion{0.0};
  std::shared_ptr<GeneFeature> gene_ptr;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct VepInfo {

  size_t male_lof{0};          // Loss of gene function in the (B) chromosome.
  size_t female_lof{0};        // Los of function in the female (A) chromosome.
  size_t male_high_effect{0};
  size_t female_high_effect{0};
  size_t male_moderate_effect{0};
  size_t female_moderate_effect{0};

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ClinvarInfo {

  std::string rsid;          // The variant rs identifier from dbSNP
  std::string clnsig;        // Textual clinical significance.
  std::string clndn;         // Clinical description.
  std::string clnisdb;       // Database information about the condition
  std::shared_ptr<const Variant> variant_ptr;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GenomeMutation {

public:

  GenomeMutation() = default;
  ~GenomeMutation() = default;

  // This analysis is performed first
  bool genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_reference);
  // Then this analysis.
  bool variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr,
                       const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                       const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                       const std::shared_ptr<const GenomePEDData>& ped_data);
  // Output to file.
  bool writeOutput(const std::string& out_file, char output_delimiter) const;

private:

  std::vector<GeneMutation> gene_vector_;

  // Vep fields.
  constexpr static const char* LOF_VEP_FIELD = "LoF";
  constexpr static const char* LOF_HC_VALUE = "HC";
  constexpr static const char* IMPACT_VEP_FIELD = "IMPACT";
  constexpr static const char* IMPACT_MODERATE_VALUE = "MODERATE";
  constexpr static const char* IMPACT_HIGH_VALUE = "HIGH";

  // Clinvar fields.
  constexpr static const char* CLINVAR_RS_FIELD = "RS";
  constexpr static const char* CLINVAR_CLNDN_FIELD = "CLNDN";
  constexpr static const char* CLINVAR_CLNSIG_FIELD = "CLNSIG";
  constexpr static const char* CLINVAR_CLNDISDB_FIELD = "CLNDISDB";
  constexpr static const char* CLINVAR_PATH_SIGNIF = "PATH";
  constexpr static const char* CLINVAR_RISK_SIGNIF = "RISK";


  void writeHeader(std::ostream& out_file, char output_delimiter) const;

  bool pedAnalysis( GeneMutation& gene_mutation,
                    const GenomeId_t& genome_id,
                    size_t data_count,
                    const std::shared_ptr<const GenomePEDData>& ped_data);

  std::shared_ptr<const ContigDB> getGeneContig( const std::shared_ptr<const ContigDB>& contig_ptr,
                                                 const GeneMutation& gene_mutation);

  std::shared_ptr<const ContigDB> getGeneSpan(const std::shared_ptr<const ContigDB>& contig_ptr,
                                              const GeneMutation& gene_mutation);

  std::shared_ptr<const ContigDB> getGeneExon(const std::shared_ptr<const ContigDB>& contig_ptr,
                                              const GeneMutation& gene_mutation);

  GeneMutation geneSpanAnalysis( const std::shared_ptr<const PopulationDB>& population_ptr,
                                 const std::shared_ptr<const PopulationDB>& unphased_population_ptr,
                                 const std::shared_ptr<const PopulationDB>& clinvar_population_ptr,
                                 const std::shared_ptr<const GenomePEDData>& ped_data,
                                 GeneMutation gene_mutation);

  VepInfo geneSpanVep( const std::shared_ptr<const ContigDB>& span_contig,
                       const std::shared_ptr<const PopulationDB>& unphased_population_ptr);

  size_t VepCount( const std::shared_ptr<const ContigDB>& vep_contig,
                   const std::string& vep_field_ident,
                   const std::string& vep_field_value);

  std::shared_ptr<const ContigDB> getClinvar( const std::shared_ptr<const ContigDB>& span_contig,
                                              const std::shared_ptr<const PopulationDB>& unphased_population_ptr);

  static std::vector<ClinvarInfo> clinvarInfo(const std::shared_ptr<const ContigDB>& clinvar_contig_ptr);

  static std::vector<ClinvarInfo> filterPathClinvar(const std::vector<ClinvarInfo>& clinvar_vector);

  static std::string clinvarConcatDesc(const std::vector<ClinvarInfo>& clinvar_vector);


};



} // namespace




#endif //KGL_KGL_ANALYSIS_MUTATION_GENE_H
