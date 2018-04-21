//
// Created by kellerberrin on 10/11/17.
//

#include "kgl_phylogenetic_env.h"
#include "kgl_phylogenetic_app.h"
#include "kgl_variant_factory.h"
#include "kgl_phylogenetic_app_analysis.h"


namespace kgl = kellerberrin::genome;


kgl::PhylogeneticApp::PhylogeneticApp(const kgl::Phylogenetic& args) {

  // Create a genome database object.
  std::shared_ptr<const GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                             args.gffFile,
                                                                                             args.gafFile,
                                                                                             args.aminoTranslationTable);
  std::shared_ptr<VCFPopulation> vcf_population_ptr(std::make_shared<VCFPopulation>());

  // For all VCF files, read in the variants.
  for (const auto& file : args.fileList) {

    kgl::VariantFactory().readVCFVariants(genome_db_ptr,
                                          vcf_population_ptr,
                                          file.genome_name,
                                          file.file_name,
                                          args.readQuality,
                                          args.variantQuality,
                                          args.minCount,
                                          args.minProportion);

  }

  // Analyze the data.
  kgl::PhylogeneticAnalysis::performAnalysis(args, genome_db_ptr, vcf_population_ptr);

}


