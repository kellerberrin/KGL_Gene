//
// Created by kellerberrin on 17/10/17.
//

#include "kgl_minority_env.h"
#include "kgl_genome_db.h"
#include "kgl_gff_fasta.h"
#include "kgl_process_sam.h"
#include "kgl_genome_analysis.h"
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::GenomeVariant> getSNPVariants(kgl::Logger& log,
                                                   std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr,
                                                   const std::string& file_name,
                                                   unsigned char read_quality,
                                                   long min_count,
                                                   double min_proportion) {

  // Create a data block to hold the read data.
  std::shared_ptr<kgl::ContigCountData> count_data_ptr(std::make_shared<kgl::ContigCountData>());

  // Register with the genome database to setup the contig data blocks.
  genome_db_ptr->registerContigData(count_data_ptr);

  // Attach a SAM reader to the contig data block and read in the SAM file.
  kgl::SamCountReader(count_data_ptr, log).readSAMFile(file_name, read_quality);

  // Generate simple SNPs.
  std::shared_ptr<kgl::GenomeVariant> variant_ptr
  = kgl::GenomeAnalysis().simpleSNPVariants<kgl::NucleotideColumn_DNA5>(count_data_ptr, genome_db_ptr);

  // Filter for read count.
  variant_ptr = variant_ptr->filterVariants(kgl::ReadCountFilter(min_count));
  // Filter for read proportion.
  variant_ptr = variant_ptr->filterVariants(kgl::MutantProportionFilter(min_proportion));

  return variant_ptr;

}


kgl::MinorityExecEnv::Application::Application(kgl::Logger& log, const kgl::MinorityArgs& args) {

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = kgl::ParseGffFasta(log).readFastaGffFile(args.fastaFile,
                                                                                                args.gffFile);
  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  // Generate mutant filtered simple SNPs.
  std::shared_ptr<kgl::GenomeVariant> mutant_variant_ptr = getSNPVariants(log,
                                                                          genome_db_ptr,
                                                                          args.mutantFile,
                                                                          args.readQuality,
                                                                          args.mutantMinCount,
                                                                          args.mutantMinProportion);
  std::cout << *mutant_variant_ptr;

  if (args.parentFile.length() > 0) {

    // Generate parent filtered simple SNPs.
    std::shared_ptr<kgl::GenomeVariant> parent_variant_ptr = getSNPVariants(log,
                                                                            genome_db_ptr,
                                                                            args.parentFile,
                                                                            args.readQuality,
                                                                            args.parentMinCount,
                                                                            args.parentMinProportion);
    std::cout << *parent_variant_ptr;


    // Union between variants in the mutant and in the parent.
    std::shared_ptr<kgl::GenomeVariant> intersect_variant_ptr = mutant_variant_ptr->Intersection(parent_variant_ptr);

    std::cout << *intersect_variant_ptr;

    // Union between variants in the mutant and in the parent.
    std::shared_ptr<kgl::GenomeVariant> union_variant_ptr = mutant_variant_ptr->Union(parent_variant_ptr);

    std::cout << *union_variant_ptr;

    // Subtract any variants from the union that are in the parent.
    std::shared_ptr<kgl::GenomeVariant> diff_variant_ptr = union_variant_ptr->Difference(parent_variant_ptr);

    std::cout << *diff_variant_ptr;

  }

}




