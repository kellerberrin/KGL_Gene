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

  // Command line arguments
  const Phylogenetic& args = getArgs();
  // XML program run options
  const RuntimeProperties& runtime_options = getRuntimeOptions();

  std::shared_ptr<kgl::GenomeCollection> genome_collection = GenomeCollection::createGenomeCollection(runtime_options);

  // A map of all active genomes.
  GenomeMap genome_map;

  // Define the organism as the pf3D7 strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> Pf3D7_aux_file_list;
  // Create the 3D7 database.
  std::shared_ptr<GenomeDatabase> genome_3D7_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::Pf3D7_41_, Pf3D7_aux_file_list);
  genome_map[genome_3D7_ptr->genomeId()] = genome_3D7_ptr;

  // Define the organism as the pfHB3 strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> PfHB3_aux_file_list;
  // Create the 3D7 database.
  std::shared_ptr<GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::PfHB3_41_, PfHB3_aux_file_list);
  genome_map[genome_db_ptr->genomeId()] = genome_db_ptr;

  // Define the organism as the pfIT strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> PfIT_aux_file_list;
  // Create the 3D7 database.
  genome_db_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::PfIT_41_, PfIT_aux_file_list);
  genome_map[genome_db_ptr->genomeId()] = genome_db_ptr;

  // Define the organism as the pfGB4 strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> PfGB4_aux_file_list;
  // Create the 3D7 database.
  genome_db_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::PfGB4_41_, PfGB4_aux_file_list);
  genome_map[genome_db_ptr->genomeId()] = genome_db_ptr;

  // Define the organism as the pfDd2 strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> PfDd2_aux_file_list;
  // Create the 3D7 database.
  genome_db_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::PfDd2_41_, PfDd2_aux_file_list);
  genome_map[genome_db_ptr->genomeId()] = genome_db_ptr;

  // Define the organism as the pf7G8 strain, PlasmoDB version 41.
  // Specify any auxiliary genome database files.
  const std::vector<std::string> Pf7G8_aux_file_list;
  // Create the 3D7 database.
  genome_db_ptr = GenomeDatabase::createGenomeDatabase(runtime_options, kgl::PhylogeneticAnalysis::Pf7G8_41_, Pf7G8_aux_file_list);
  genome_map[genome_db_ptr->genomeId()] = genome_db_ptr;

  // Create a phased population object.
  std::shared_ptr<PhasedPopulation> population_ptr(std::make_shared<PhasedPopulation>(kgl::PhylogeneticAnalysis::Pf3D7_41_));

  // Create an unphased population object.
  std::shared_ptr<UnphasedPopulation> unphased_population_ptr(std::make_shared<UnphasedPopulation>());

  // Write Filtered Unphased Heterozygous Statistics
  HeterozygousStatistics heterozygous_statistics;

  std::vector<std::string> vcf_list;
  runtime_options.getVCFFiles(vcf_list);



  // For all VCF files, read in the variants.
  for (auto vcf_file : vcf_list) {

    // Clear the unphased population object.
    unphased_population_ptr->clear();

    // Read variants.
    VariantFactory().readVCFVariants(genome_3D7_ptr, unphased_population_ptr, vcf_file);

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
      GenomePhasing::fileHaploidPhasing(mixture_file, ploidy, filtered_unphased_ptr, genome_3D7_ptr, population_ptr);

    } else {

      // No mixture file, so assume all genomes are unmixed and clonal
      GenomePhasing::haploidPhasing(ploidy, filtered_unphased_ptr, genome_3D7_ptr, population_ptr);

    }

  }

  // Write the unphased hetero/homozygous statistics.
  std::string heterozygous_file = kgl::Utility::filePath("HeterozygousAnalysis", args.workDirectory) + ".csv";
  heterozygous_statistics.writeHeterozygousStatistics(heterozygous_file, ',');

  // Analyze the data.
  kgl::PhylogeneticAnalysis(runtime_options, genome_map, unphased_population_ptr, population_ptr).performAnalysis(args.analysisType);

}

