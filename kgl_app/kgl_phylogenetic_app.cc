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


namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant>
kgl::PhylogeneticApp::getGenomeVariants(std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                        const std::string& file_name,
                                        bool vcf_is_gatk,
                                        const std::string& genome_name,
                                        kgl::Phred_t read_quality,
                                        kgl::Phred_t variant_quality,
                                        long min_count,
                                        double min_proportion) {

  // Read in the SAM/VCF file variants
  std::shared_ptr<const kgl::GenomeVariant> all_variant_ptr = kgl::VariantFactory().createVariants(genome_db_ptr,
                                                                                                   genome_name,
                                                                                                   file_name,
                                                                                                   vcf_is_gatk,
                                                                                                   read_quality,
                                                                                                   variant_quality,
                                                                                                   min_count,
                                                                                                   min_proportion);

  // Filter on  quality >= 5.
  read_quality = read_quality < 5 ? 5 : read_quality;

  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = all_variant_ptr->filterVariants(kgl::QualityFilter(read_quality));

  kgl::ExecEnv::log().info("Filtered for quality: {}, Genome: {} has: {} variants", read_quality, genome_name, filter_ptr->size());

  // Return the genome variants.
  return filter_ptr;

}


kgl::PhylogeneticApp::PhylogeneticApp(const kgl::Phylogenetic& args) {


  // Create a genome database object.
  std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr = GenomeDatabase::createGenomeDatabase(args.fastaFile,
                                                                                                  args.gffFile,
                                                                                                  args.gafFile,
                                                                                                  args.aminoTranslationTable);

  // Create a population object.
  std::shared_ptr<kgl::PopulationVariant> pop_variant_ptr(std::make_shared<kgl::PopulationVariant>("Falciparum"));

  // For all organisms
  for (const auto& file : args.fileList) {

    // Generate organism genome variants.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getGenomeVariants(genome_db_ptr,
                                                                              file.file_name,
                                                                              args.vcfAllGATK,
                                                                              file.genome_name,
                                                                              args.readQuality,
                                                                              args.variantQuality,
                                                                              args.minCount,
                                                                              args.minProportion);

    // Store the organism variants in the population object.
    pop_variant_ptr->addGenomeVariant(variant_ptr);



  }

  // Analyze the population.
  performAnalysis(args, genome_db_ptr, pop_variant_ptr);

}


