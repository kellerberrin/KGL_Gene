//
// Created by kellerberrin on 17/10/17.
//
#include <sstream>

#include "kgl_minority_env.h"
#include "kgl_genome_db.h"
#include "kgl_gff_fasta.h"
#include "kgl_sam_process.h"
#include "kgl_variant_evidence.h"
#include "kgl_filter.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant> getSNPVariants(kgl::Logger& log,
                                                         std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                         const std::string& file_name,
                                                         unsigned char read_quality,
                                                         long min_count,
                                                         double min_proportion) {


  // Read in the SAM file.
  std::shared_ptr<const kgl::ContigCountData> count_data_ptr = kgl::SamCountReader(log).readSAMFile(genome_db_ptr,
                                                                                                    file_name,
                                                                                                    read_quality);
  // Generate SNP variants.
  std::shared_ptr<const kgl::GenomeVariant> variant_ptr = kgl::VariantAnalysis().SNPVariants(count_data_ptr,
                                                                                             genome_db_ptr);
  // Generate contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> codon_delete_ptr = kgl::VariantAnalysis().codonDelete(variant_ptr,
                                                                                                  count_data_ptr,
                                                                                                  genome_db_ptr);
  variant_ptr = variant_ptr->filterVariants(kgl::NotFilter<kgl::DeleteSNPFilter>());
  // Filter for CDS membership.
  variant_ptr = variant_ptr->filterVariants(kgl::InCDSFilter(genome_db_ptr));
  std::cout << *variant_ptr;
  // Filter for read count.
  variant_ptr = variant_ptr->filterVariants(kgl::ReadCountFilter(min_count));
  // Filter for read proportion.
  variant_ptr = variant_ptr->filterVariants(kgl::MutantProportionFilter(min_proportion));
  // Add in compound deletes.
  variant_ptr = variant_ptr->Union(codon_delete_ptr);
  // Filter for PfATP4
  variant_ptr = variant_ptr->filterVariants(kgl::GeneFilter("PF3D7_1211900",genome_db_ptr));

  return variant_ptr;

}


kgl::MinorityExecEnv::Application::Application(kgl::Logger& log, const kgl::MinorityArgs& args) {

  log.SetVerbose(args.verbose);

  // Create a genome database object.
  std::shared_ptr<kgl::GenomeDatabase> genome_db_ptr = kgl::ParseGffFasta(log).readFastaGffFile(args.fastaFile,
                                                                                                args.gffFile);
  // Set the amino translation table
  genome_db_ptr->setTranslationTable(args.aminoTranslationTable);

  // Wire-up the genome database.
  genome_db_ptr->createVerifyGenomeDatabase();

  std::vector<std::shared_ptr<const kgl::GenomeVariant>> variant_vector;
  if (not args.fileList.empty()) {

    for(auto sam_file : args.fileList) {

      // Generate mutant filtered simple SNPs.
      std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getSNPVariants(log,
                                                                             genome_db_ptr,
                                                                             sam_file,
                                                                             args.readQuality,
                                                                             args.mutantMinCount,
                                                                             args.mutantMinProportion);
      variant_vector.push_back(variant_ptr);
      std::ostringstream ss;
      ss << *variant_ptr;
      ExecEnv::log().info("Genome: {}\nPF3D7_1211900 variants\n{}", variant_ptr->genomeId(), ss.str());

    }

  }


  if (args.mutantFile.length() > 0) {
    // Generate mutant filtered simple SNPs.
    std::shared_ptr<const kgl::GenomeVariant> mutant_variant_ptr = getSNPVariants(log,
                                                                                  genome_db_ptr,
                                                                                  args.mutantFile,
                                                                                  args.readQuality,
                                                                                  args.mutantMinCount,
                                                                                  args.mutantMinProportion);
//    std::cout << *mutant_variant_ptr;

  }

  if (args.parentFile.length() > 0) {

    // Generate parent filtered simple SNPs.
    std::shared_ptr<const kgl::GenomeVariant> parent_variant_ptr = getSNPVariants(log,
                                                                                  genome_db_ptr,
                                                                                  args.parentFile,
                                                                                  args.readQuality,
                                                                                  args.parentMinCount,
                                                                                  args.parentMinProportion);
//    std::cout << *parent_variant_ptr;


    // Union between variants in the mutant and in the parent.
//    std::shared_ptr<kgl::GenomeVariant> intersect_variant_ptr = mutant_variant_ptr->Intersection(parent_variant_ptr);

//    std::cout << *intersect_variant_ptr;

    // Union between variants in the mutant and in the parent.
//    std::shared_ptr<kgl::GenomeVariant> union_variant_ptr = mutant_variant_ptr->Union(parent_variant_ptr);

//    std::cout << *union_variant_ptr;

    // Subtract any variants from the union that are in the parent.
//    std::shared_ptr<kgl::GenomeVariant> diff_variant_ptr = union_variant_ptr->Difference(parent_variant_ptr);

//    std::cout << *diff_variant_ptr;

  }

}




