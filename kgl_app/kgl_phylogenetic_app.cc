//
// Created by kellerberrin on 10/11/17.
//

#include "../kgl_mixture/kgl_variant_phase.h"
#include "kgl_phylogenetic_app.h"
#include "kgl_variant_factory.h"
#include "kgl_phylogenetic_app_analysis.h"
#include "kgl_unphased_analysis.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;


void kgl::PhylogeneticExecEnv::executeApp() {

  const Phylogenetic& args = getArgs();
  const RuntimeOptions& runtime_options = getRuntimeOptions();

  // Create a phased population object.
  std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>("Falciparum"));

  // Create a genome database object.
  std::shared_ptr<const GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                             args.gffFile,
                                                                                             args.gafFile,
                                                                                             args.aminoTranslationTable);

  HeterozygousStatistics heterozygous_statistics;
  // For all VCF files, read in the variants.
  for (const auto& file : args.fileList) {

    std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>());

    kgl::VariantFactory().readVCFVariants(genome_db_ptr,
                                          unphased_population_ptr,
                                          file.genome_name,
                                          file.file_name,
                                          args.variantQuality,
                                          args.minCount,
                                          args.minProportion);

    // Basic statistics to output
    unphased_population_ptr->popStatistics();
    // Filtered Unphased Heterozygous Statistics
    std::shared_ptr<UnphasedPopulation> filtered_unphased_ptr = unphased_population_ptr->filterVariants(CountFilter(30));
    // If the mixture file is defined then read it and phase the variants.
    std::string mixture_file;
    if (runtime_options.getMixtureFile(mixture_file)) {

      GenomePhasing::fileHaploidPhasing(mixture_file, 2, unphased_population_ptr, genome_db_ptr, population_ptr);

    }
    // Analyze the data.
    kgl::PhylogeneticAnalysis::performAnalysis(args, runtime_options, genome_db_ptr, filtered_unphased_ptr, population_ptr);
    // Process Filtered Unphased Heterozygous Statistics
    heterozygous_statistics.heterozygousStatistics(filtered_unphased_ptr);

  }

  // Write Filtered Unphased Heterozygous Statistics
  std::string heterozygous_file = kgl::Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
  heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',');

}


