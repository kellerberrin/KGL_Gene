//
// Created by kellerberrin on 31/1/21.
//

#ifndef KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
#define KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H

#include "kgl_variant_db_population.h"
#include "kgl_analysis_mutation_gene_ethnic.h"



namespace kellerberrin::genome {   //  organization::project level namespace


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantPhaseStats {

public:

  VariantPhaseStats() = default;
  VariantPhaseStats(const VariantPhaseStats &) = default;
  ~VariantPhaseStats() = default;

  VariantPhaseStats &operator=(const VariantPhaseStats &) = default;

  bool phaseStatistics(const std::shared_ptr<const ContigDB>& contig_ptr);
  static void writeHeader(std::ostream& out_file, char output_delimiter);


    // Variant phase
  size_t phase_male_{0};   // Allele(s) is present on the male_ phase
  size_t phase_female_{0};  // Allele(s) is present on the female_ phase
  size_t phase_hom_{0};  // Allele(s) is present on both phases. Male and female_.
  size_t phase_either_{0};  // Allele(s) is present on either phase. Male or female_.

private:


};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneClinvar {

public:

  GeneClinvar() = default;

  GeneClinvar(const GeneClinvar &) = default;

  ~GeneClinvar() = default;

  GeneClinvar &operator=(const GeneClinvar &) = default;


  // Variant phase
  VariantPhaseStats phase;
  // Text description of problems with this gene.
  std::set<std::string> clinvar_desc;
  // Superpopulation and sex breakdown.
  GeneEthnicitySex results;

};


} // namespace


#endif // KGL_ANALYSIS_MUTATION_GENE_CLINVAR_H
