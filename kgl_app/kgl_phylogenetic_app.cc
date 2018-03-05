//
// Created by kellerberrin on 10/11/17.
//

#include <sstream>
#include "kgl_utility.h"
#include "kgl_genome_types.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_phylogenetic_app.h"
#include "kgl_variant_factory.h"
#include "kgl_filter.h"
#include "kgl_ploidy_analysis.h"


namespace kgl = kellerberrin::genome;


kgl::PhylogeneticApp::PhylogeneticApp(const kgl::Phylogenetic& args) {

  // Create a genome database object.
  std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                                  args.gffFile,
                                                                                                  args.gafFile,
                                                                                                  args.aminoTranslationTable);
  // Create a population object.
  std::shared_ptr<kgl::PopulationVariant> pop_variant_ptr(std::make_shared<kgl::PopulationVariant>("Falciparum"));

  std::shared_ptr<PloidyAnalysis> ploidy_ptr(std::make_shared<PloidyAnalysis>("PloidyAnalysis"));

  // For all organisms
  for (const auto& file : args.fileList) {

/*
    kgl::VariantFactory().createVariants(genome_db_ptr,
                                         pop_variant_ptr,
                                         file.genome_name,
                                         file.file_name,
                                         args.readQuality,
                                         args.variantQuality,
                                         args.minCount,
                                         args.minProportion);
*/

    kgl::VariantFactory().createVariants(genome_db_ptr,
                                         ploidy_ptr,
                                         file.genome_name,
                                         file.file_name,
                                         args.readQuality,
                                         args.variantQuality,
                                         args.minCount,
                                         args.minProportion);

  }

  std::string ploidy_file = kgl::Utility::filePath("PloidyAnalysis", args.workDirectory) + ".csv";
  ploidy_ptr->writePloidyResults(ploidy_file, PloidyAnalysis::CSV_DELIMITER_);

  // Analyze the population.
  performAnalysis(args, genome_db_ptr, pop_variant_ptr);

}


