//
// Created by kellerberrin on 10/11/17.
//

#include <sstream>
#include "kgl_utility.h"
#include "kgl_library/kgl_genome_types.h"
#include "kgl_phylogenetic_env.h"
#include "kgl_library/kgl_genome_db.h"
#include "kgl_library/kgl_gff_fasta.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory.h"
#include "kgl_library/kgl_variant_factory_compound.h"
#include "kgl_library/kgl_filter.h"
#include "kgl_phylogenetic_analysis.h"
#include "kgl_statistics.h"
#include "kgl_phylogenetic_gene.h"

namespace kgl = kellerberrin::genome;


// pfATP4 drug target ATP4 sodium pump.
#define PFATP4_MINORITY_CONTIG "chr12"
#define PFATP4_MINORITY_GENE "PF3D7_1211900"
#define PFATP4_MINORITY_SEQUENCE "rna_PF3D7_1211900-1"

#define PFATP4_MALAWI_CONTIG "Pf3D7_12_v3"
#define PFATP4_MALAWI_GENE "PF3D7_1211900"
#define PFATP4_MALAWI_SEQUENCE "PF3D7_1211900.1"

// Erythrocyte membrane (Blood Cell Surface) protein (busy Malawi)
#define _BCS_CONTIG "Pf3D7_12_v3"
#define _BCS_GENE "PF3D7_1255200"
#define  BCS_SEQUENCE "PF3D7_1255200.1"

// Mutant rich Rifin (Malawi)
#define RIFIN_3_CONTIG "Pf3D7_01_v3"
#define RIFIN_3_GENE "PF3D7_0101900"
#define RIFIN_3_SEQUENCE "PF3D7_0101900.1"

// Mutant rich Rifin (Malawi)
#define RIFIN_4_CONTIG "Pf3D7_08_v3"
#define RIFIN_4_GENE "PF3D7_0808900"
#define RIFIN_4_SEQUENCE "PF3D7_0808900.1"

// Rifin very similar mutations (Malawi).
#define RIFIN_1_CONTIG "Pf3D7_07_v3"
#define RIFIN_1_GENE "PF3D7_0711700"
#define RIFIN_1_SEQUENCE "PF3D7_0711700.1"

#define RIFIN_2_CONTIG "Pf3D7_07_v3"
#define RIFIN_2_GENE  "PF3D7_0712600"
#define RIFIN_2_SEQUENCE "PF3D7_0712600.1"

// S-Antigen very busy 5-prime and 3-prime regions (Malawi).
#define S_ANTIGEN_CONTIG "Pf3D7_10_v3"
#define S_ANTIGEN_GENE "PF3D7_1035200"
#define S_ANTIGEN_SEQUENCE "PF3D7_1035200.1"

#define ACTIVE_CONTIG PFATP4_MALAWI_CONTIG
#define ACTIVE_GENE PFATP4_MALAWI_GENE
#define ACTIVE_SEQUENCE PFATP4_MALAWI_SEQUENCE





std::shared_ptr<const kgl::GenomeVariant> getGenomeVariants(std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                            const std::string& file_name,
                                                            const std::string& genome_name,
                                                            kgl::Phred_t read_quality,
                                                            kgl::Phred_t variant_quality,
                                                            long min_count,
                                                            double min_proportion,
                                                            const std::string& workDirectory) {

  // Read in the SAM file variants
  std::shared_ptr<const kgl::GenomeVariant> all_variant_ptr = kgl::VariantFactory().createVariants(genome_db_ptr,
                                                                                                   genome_name,
                                                                                                   file_name,
                                                                                                   read_quality,
                                                                                                   variant_quality,
                                                                                                   min_count,
                                                                                                   min_proportion);


  // Filter on sequence and quality >= 5.
  std::shared_ptr<const kgl::GenomeVariant> filter_ptr = all_variant_ptr->filterVariants(kgl::QualityFilter(5));
  // Filter on contig

  kgl::ExecEnv::log().info("Filtered for quality: {}, Genome: {} has: {} variants", 5, genome_name, filter_ptr->size());

  // Return the genome variants.
  return filter_ptr;

}


kgl::PhylogeneticExecEnv::Application::Application(kgl::Logger& log, const kgl::Phylogenetic& args) {

  log.SetVerbose(args.verbose);

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = kgl::ParseGffFasta().readFastaGffFile(args.fastaFile,
                                                                                             args.gffFile);
  // Set the amino translation table
  genome_db_ptr->setTranslationTable(args.aminoTranslationTable);

  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  // Create a population variant object.
  std::shared_ptr<kgl::PopulationVariant> pop_variant_ptr(std::make_shared<kgl::PopulationVariant>("Malawi-PRJNA173723"));

  // For all organisms
  for (const auto& file : args.fileList) {

    // Generate all genome variants.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getGenomeVariants(genome_db_ptr,
                                                                              file.file_name,
                                                                              file.genome_name,
                                                                              args.readQuality,
                                                                              args.variantQuality,
                                                                              args.minCount,
                                                                              args.minProportion,
                                                                              args.workDirectory);

    std::shared_ptr<GenomeVariant> filter_ptr = variant_ptr->filterVariants(kgl::GeneFilter("PF3D7_0900600"));
    std::cout << filter_ptr;

    // Store the genome variant pointer
    pop_variant_ptr->addGenomeVariant(variant_ptr);



  }

//  ApplicationAnalysis::outputSequenceCSV(args.outCSVFile, genome_db_ptr, pop_variant_ptr);

  GeneAnalysis::mutateGenomeRegion("vcf_SRR609073", "Pf3D7_09_v3", 46890, 30, pop_variant_ptr, genome_db_ptr);

  // Perform population analysis
  // (disabled to save CPU time)
  std::string newick_file = Utility::filePath("file.genome_name+ UPGMA_newick", args.workDirectory) + ".txt";
//  kgl::PhylogeneticAnalysis::UPGMA(newick_file, population_stats_ptr);

}


