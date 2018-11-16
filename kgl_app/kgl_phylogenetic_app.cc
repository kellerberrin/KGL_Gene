//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_gff_fasta.h"
#include "kgl_variant_phase.h"
#include "kgl_phylogenetic_app.h"
#include "kgl_variant_factory.h"
#include "kgl_phylogenetic_app_analysis.h"
#include "kgl_distribution_analysis.h"
#include "kgl_filter.h"
#include "kgl_properties.h"


namespace kgl = kellerberrin::genome;


void kgl::PhylogeneticExecEnv::executeApp() {

  const Phylogenetic& args = getArgs();
  const RuntimeProperties& runtime_options = getRuntimeOptions();


  // Create a phased population object.
  std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>("Falciparum"));

  // Create an unphased population object.
  std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>());

  // Get the genome database runtime parameters.
  std::string fasta_file, gff_file, gaf_file, tss_file, amino_translation_table;
  runtime_options.getGenomeDBFiles(fasta_file, gff_file, gaf_file, tss_file, amino_translation_table);

  // Create a genome database object.
  std::shared_ptr<const GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(fasta_file,
                                                                                             gff_file,
                                                                                             tss_file,
                                                                                             gaf_file,
                                                                                             amino_translation_table);

  // Write Filtered Unphased Heterozygous Statistics
  HeterozygousStatistics heterozygous_statistics;

  std::vector<std::string> vcf_list;
  runtime_options.getVCFFiles(vcf_list);

  // For all VCF files, read in the variants.
  for (const auto& file : vcf_list) {

    // Clear the unphased population object.
    unphased_population_ptr->clear();

    // Read variants.
    VariantFactory().readVCFVariants(genome_db_ptr, unphased_population_ptr, file);

    // Basic statistics to output
    // unphased_population_ptr->popStatistics();

    // Filter unphased variants for minimum read statistics.
    std::shared_ptr<UnphasedPopulation> filtered_unphased_ptr = unphased_population_ptr->filterVariants(AndFilter(DPCountFilter(30), RefAltCountFilter(30)));

    // Process Filtered Unphased Heterozygous Statistics
    heterozygous_statistics.heterozygousStatistics(filtered_unphased_ptr);

    // Get the VCF ploidy (need not be the organism ploidy).
    size_t ploidy = runtime_options.getVCFPloidy();

    // If the mixture file is defined and exists then read it and phase the variants.
    std::string mixture_file;
    if (runtime_options.getMixtureFile(mixture_file)) {

      // The mixture file indicates which VCF samples have Complexity Of Infection (COI), only clonal infections are used.
      GenomePhasing::fileHaploidPhasing(mixture_file, ploidy, filtered_unphased_ptr, genome_db_ptr, population_ptr);

    } else {

      // No mixture file, so assume all genomes are unmixed and clonal
      GenomePhasing::haploidPhasing(ploidy, filtered_unphased_ptr, genome_db_ptr, population_ptr);

    }

  }

  // Write the unphased hetero/homozygous statistics.
  std::string heterozygous_file = kgl::Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
  heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',');

  // Analyze the data.
  kgl::PhylogeneticAnalysis(runtime_options, genome_db_ptr, unphased_population_ptr, population_ptr).performAnalysis(args.analysisType);

}


