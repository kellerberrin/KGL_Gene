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
  const PropertyTree& runtime_options = getRuntimeOptions();


  // Create a phased population object.
  std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>("Falciparum"));

  // Create an unphased population object.
  std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>());

  std::string tss_file;
  runtime_options.getProperty("tssFile", tss_file);
  tss_file = kgl::Utility::filePath(tss_file, args.workDirectory);

  // Create a genome database object.
  std::shared_ptr<const GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                             args.gffFile,
                                                                                             tss_file,
                                                                                             args.gafFile,
                                                                                             args.aminoTranslationTable);

  // Write Filtered Unphased Heterozygous Statistics
  HeterozygousStatistics heterozygous_statistics;

  // For all VCF files, read in the variants.
  for (const auto& file : args.fileList) {

    // clear the unphased population object.
    unphased_population_ptr->clear();

    // Read variants.
    VariantFactory().readVCFVariants(genome_db_ptr, unphased_population_ptr, file.file_name);

    // Basic statistics to output
    // unphased_population_ptr->popStatistics();

    // Filtered Unphased Heterozygous Statistics
    std::shared_ptr<UnphasedPopulation> filtered_unphased_ptr = unphased_population_ptr->filterVariants(AndFilter(DPCountFilter(30), RefAltCountFilter(30)));

    // Process Filtered Unphased Heterozygous Statistics
    heterozygous_statistics.heterozygousStatistics(filtered_unphased_ptr);

    // If the mixture file is defined then read it and phase the variants.
    std::string mixture_file;
    if (runtime_options.checkProperty("mixtureFile")) {

      runtime_options.getProperty("mixtureFile", mixture_file);
      GenomePhasing::fileHaploidPhasing(mixture_file, 2, filtered_unphased_ptr, genome_db_ptr, population_ptr);

    } else {

      // No mixture file, so assume all genomes are unmixed.
      GenomePhasing::haploidPhasing(2, filtered_unphased_ptr, genome_db_ptr, population_ptr);

    }

  }

  // Write the unphased hetero/homozygous statistics.
  std::string heterozygous_file = kgl::Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
  heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',');

  // Analyze the data.
  kgl::PhylogeneticAnalysis::performAnalysis(args, runtime_options, genome_db_ptr, unphased_population_ptr, population_ptr);

}


