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
  size_t variant_count{0};
  size_t male_phase{0};
  size_t female_phase{0};
  size_t EAS_variant_count{0};
  size_t EUR_variant_count{0};
  size_t AFR_variant_count{0};
  size_t AMR_variant_count{0};
  size_t SAS_variant_count{0};
  size_t genome_count{0};   // Total number of genomes.
  size_t male_variant{0};    // Males that have variants for this gene.
  size_t female_variant{0};  // Females that have variants for thgis gene.
  size_t genome_variant{0};  // Number of genomes that contain variants for this gene.
  size_t homozygous{0};
  size_t heterozygous{0};
  double indel{0.0};
  double transition{0.0};
  double transversion{0.0};

};


class GenomeMutation {

public:

  GenomeMutation() = default;
  ~GenomeMutation() = default;

  // This analysis is performed first
  bool genomeAnalysis( const std::shared_ptr<const GenomeReference>& genome_reference);
  // Then this analysis for Homosapien.
  bool variantAnalysis100(const std::shared_ptr<const PopulationDB>& population_ptr,
                          const std::shared_ptr<const GenomePEDData>& ped_data);

  // Then this analysis for P. Falciparum.
  bool variantAnalysis(const std::shared_ptr<const PopulationDB>& population_ptr);

  [[nodiscard]] const std::vector<GeneMutation>& geneVector() const { return gene_vector_; }
  [[nodiscard]] std::vector<GeneMutation>& geneVector() { return gene_vector_; }

  bool writeOutput(const std::string& out_file, char output_delimiter) const;
  bool writeOutput100(const std::string& out_file, char output_delimiter) const;

private:

  std::vector<GeneMutation> gene_vector_;

  void writeHeader(std::ostream& out_file, char output_delimiter) const;
  void writeHeader100(std::ostream& out_file, char output_delimiter) const;

};



} // namespace




#endif //KGL_KGL_ANALYSIS_MUTATION_GENE_H
