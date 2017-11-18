//
// Created by kellerberrin on 10/11/17.
//

#include <sstream>
#include "kgl_phylogenetic_env.h"
#include "kgl_genome_db.h"
#include "kgl_gff_fasta.h"
#include "kgl_sam_process.h"
#include "kgl_variant_evidence.h"
#include "kgl_filter.h"
#include "kgl_phylogenetic_analysis.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant> getSNPVariants(kgl::Logger& log,
                                                         std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                                         const std::string& file_name,
                                                         const std::string& genome_name,
                                                         unsigned char read_quality,
                                                         long min_count,
                                                         double min_proportion) {

  // Read in the SAM file.
  std::shared_ptr<const kgl::ContigCountData> count_data_ptr = kgl::SamCountReader(log).readSAMFile(genome_db_ptr,
                                                                                                    file_name,
                                                                                                    read_quality);
  // Generate SNP variants.
  std::shared_ptr<const kgl::GenomeVariant> variant_ptr = kgl::VariantAnalysis().SNPVariants(genome_name,
                                                                                             count_data_ptr,
                                                                                             genome_db_ptr,
                                                                                             min_count,
                                                                                             min_proportion);
  // Generate contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> codon_delete_ptr = kgl::VariantAnalysis().codonDelete(variant_ptr,
                                                                                                  count_data_ptr,
                                                                                                  genome_db_ptr);
  // Disaggregated contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> disagg_ptr = codon_delete_ptr->disaggregateCompoundVariants(genome_db_ptr);
  // Remove disaggregated variants.
  variant_ptr = variant_ptr->Difference(disagg_ptr);
  // Add in contiguous deletes.
  variant_ptr = variant_ptr->Union(codon_delete_ptr);
  // Generate compound single codon variants
  std::shared_ptr<const kgl::GenomeVariant> compound_snp_ptr = kgl::VariantAnalysis().compoundSNP(variant_ptr,
                                                                                                  genome_db_ptr);
  variant_ptr = variant_ptr->filterVariants(kgl::InCDSFilter());
  std::cout << *variant_ptr;
  // Filter for PfATP4
  variant_ptr = variant_ptr->filterVariants(kgl::GeneFilter("PF3D7_1211900"));
//  variant_ptr = variant_ptr->filterVariants(kgl::SequenceFilter("rna_PF3D7_1211900-1"));
  return variant_ptr;

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

  for (const auto& file : args.fileList) {

    // Generate mutant filtered simple SNPs.
    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = getSNPVariants(log,
                                                                           genome_db_ptr,
                                                                           file.file_name,
                                                                           file.genome_name,
                                                                           args.readQuality,
                                                                           args.minCount,
                                                                           args.minProportion);


    variant_ptr->outputCSV(args.outCSVFile, VariantOutputIndex::START_1_BASED);
    std::string fasta_file_name = ExecEnv::filePath(file.genome_name, args.workDirectory) + ".fasta";
    std::string sequence_name = "PF3D7_1211900_" + file.genome_name;
    std::string read_fasta_name = ExecEnv::filePath("PF3D7_1211900_Reference", args.workDirectory) + ".fasta";


#define MUTANT_PROTEIN 1
#define MALAWI 1

#ifdef MUTANT_PROTEIN
#ifdef MALAWI

    ApplicationAnalysis::writeMutantProtein(fasta_file_name,
                                            sequence_name,
                                            "Pf3D7_12_v3",
                                            "PF3D7_1211900",
                                            "PF3D7_1211900.1",
                                            genome_db_ptr,
                                            variant_ptr);

    std::string comparison;
    if (ApplicationAnalysis::readMutantProtein(read_fasta_name,
                                           "PF3D7_1211900",
                                           "Pf3D7_12_v3",
                                           "PF3D7_1211900",
                                           "PF3D7_1211900.1",
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison)) {

      ExecEnv::log().info("PF3D7_1211900 comparison:\n{}", comparison);

    }

#else

    ApplicationAnalysis::writeMutantProtein(fasta_file_name,
                                            sequence_name,
                                            "chr12",
                                            "PF3D7_1211900",
                                            "rna_PF3D7_1211900-1",
                                            genome_db_ptr,
                                            variant_ptr);
    std::string comparison;
    if (ApplicationAnalysis::readMutantProtein(read_fasta_name,
                                           "PF3D7_1211900",
                                           "chr12",
                                           "PF3D7_1211900",
                                           "rna_PF3D7_1211900-1",
                                           genome_db_ptr,
                                           variant_ptr,
                                           comparison)) {

      ExecEnv::log().info("PF3D7_1211900 comparison:\n{}", comparison);

    }

#endif
#endif

  }


}


