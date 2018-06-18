//
// Created by kellerberrin on 27/11/17.
//

#include "kgl_utility.h"
#include "kgl_sam_process.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_bam.h"
#include "kgl_variant_single.h"
#include "kgl_variant_factory.h"
#include "kgl_variant_factory_single.h"
#include "kgl_variant_factory_compound.h"
#include "kgl_filter.h"
#include "kgl_variant_factory_vcf_impl.h"


namespace kgl = kellerberrin::genome;



void kgl::VariantFactory::readVCFVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                          std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                          const std::string& genome_name,
                                          const std::string& variant_file_name,
                                          Phred_t variant_quality,
                                          NucleotideReadCount_t min_read_count,
                                          double min_proportion) const {

  std::string file_ext = kgl::Utility::fileExtension(variant_file_name);
  std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::toupper); // convert to UC for robust comparison

  if (file_ext == VCF_FILE_EXTENSTION_) {

      ExecEnv::log().info("Processing VCF file: {}", variant_file_name);
      VcfFactory().readParseVCFVariants(vcf_population_ptr, genome_db_ptr, variant_file_name, variant_quality);

  } else {

    ExecEnv::log().error("Invalid file name: {}", variant_file_name);
    ExecEnv::log().critical("Unsupported file type: '{}' for variant calling. Must be VCF ('.vcf')", file_ext);

  }

}


void kgl::VariantFactory::readCountVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                            std::shared_ptr<PhasedPopulation> pop_variant_ptr,
                                            const std::string& genome_name,
                                            const std::string& variant_file_name,
                                            Phred_t read_quality,
                                            NucleotideReadCount_t min_read_count,
                                            double min_proportion) const {

  std::string file_ext = kgl::Utility::fileExtension(variant_file_name);
  std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::toupper); // convert to UC for robust comparison

  if (file_ext == SAM_FILE_EXTENSTION_) {

    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = createSamVariants(genome_db_ptr,
                                                                              genome_name,
                                                                              variant_file_name,
                                                                              read_quality,
                                                                              min_read_count,
                                                                              min_proportion);
    addGenome(variant_ptr, pop_variant_ptr, read_quality);

  } else if (file_ext == BAM_FILE_EXTENSTION_) {

    std::shared_ptr<const kgl::GenomeVariant> variant_ptr = createBamVariants(genome_db_ptr,
                                                                              genome_name,
                                                                              variant_file_name,
                                                                              read_quality,
                                                                              min_read_count,
                                                                              min_proportion);
    addGenome(variant_ptr, pop_variant_ptr, read_quality);

  } else {

    ExecEnv::log().error("Invalid file name: {}", variant_file_name);
    ExecEnv::log().critical("Unsupported file type: '{}' for variant calling. Must be SAM ('.sam') or BAM ('.bam')", file_ext);

  }

}



bool kgl::VariantFactory::isFileNamePrefix(const std::string& prefix, const std::string& variant_file_name) const {

  std::string file_name = kgl::Utility::fileName(variant_file_name);
  std::transform(file_name.begin(), file_name.end(), file_name.begin(), ::toupper); // convert to UC for robust comparison
  std::string file_name_prefix = file_name.substr(0, prefix.length());

  return (file_name_prefix == prefix);

}


void kgl::VariantFactory::addGenome(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                                    std::shared_ptr<PhasedPopulation> pop_variant_ptr,
                                    Phred_t read_quality) const {


// Store the organism variants in the population object.
  pop_variant_ptr->addGenomeVariant(genome_variant_ptr);

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createSamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::string& sam_file_name,
                                       Phred_t read_quality,
                                       NucleotideReadCount_t min_read_count,
                                       double min_proportion) const {

  // Read in the SAM file.
  std::shared_ptr<const kgl::ContigCountData> count_data_ptr = kgl::SamCountReader().readSAMFile(genome_db_ptr,
                                                                                                 sam_file_name,
                                                                                                 read_quality);

  ExecEnv::log().info("Processing SAM file: {}", sam_file_name);
  ExecEnv::log().info("Generating Pileup (count) variants for Genome: {}", genome_name);

  // generate raw variants.
  std::shared_ptr<const GenomeVariant> single_variant_ptr = SingleFactory().createSingleVariants(genome_name,
                                                                                                 count_data_ptr,
                                                                                                 genome_db_ptr,
                                                                                                 min_read_count,
                                                                                                 min_proportion);

  ExecEnv::log().info("Generated: {} Pileup (count) variants for Genome: {}", single_variant_ptr->variantCount(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, single_variant_ptr);

  return variant_ptr;

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createBamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::string& bam_file_name,
                                       Phred_t read_quality,
                                       NucleotideReadCount_t min_read_count,
                                       double min_proportion) const {


  ExecEnv::log().info("Processing BAM file: {}", bam_file_name);
  ExecEnv::log().info("Generating Pileup (count) variants for Genome: {}", genome_name);
  // generate snp variants.
  std::shared_ptr<const GenomeVariant> single_variant_ptr = BamFactory().readParseBam(genome_name,
                                                                                      genome_db_ptr,
                                                                                      bam_file_name,
                                                                                      read_quality,
                                                                                      min_read_count,
                                                                                      min_proportion);

  ExecEnv::log().info("Generated: {} Pileup (count) variants for Genome: {}", single_variant_ptr->variantCount(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, single_variant_ptr);

  return variant_ptr;

}


// Generate compound variants.
std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::aggregateVariants(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::shared_ptr<const kgl::GenomeVariant>& single_variant_ptr) {


  ExecEnv::log().info("Generating Compound variants for Genome: {}", genome_name);

  // Generate contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_delete_ptr = kgl::CompoundDeleteFactory().create(single_variant_ptr,
                                                                                                 genome_db_ptr);
  // Generate contiguous insertion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_insert_ptr = kgl::CompoundInsertFactory().create(single_variant_ptr,
                                                                                                 genome_db_ptr);

  // combine the compound variants
  std::shared_ptr<const kgl::GenomeVariant> compound_ptr = cmp_delete_ptr->Union(cmp_insert_ptr);

  ExecEnv::log().info("Generated: {} Compound variants for Genome: {}", compound_ptr->variantCount(), genome_name);
  ExecEnv::log().info("Combining Compound and Single variants for Genome: {}", genome_name);

  // disaggregate the compound variants.
  std::shared_ptr<const kgl::GenomeVariant> disagg_ptr = kgl::CompoundDeleteFactory().disaggregate(compound_ptr,
                                                                                                   genome_db_ptr);

  // remove the disaggregated snps from the set of raw snp variants.
  std::shared_ptr<const kgl::GenomeVariant> variant_ptr = single_variant_ptr->Difference(disagg_ptr);

  // add in the compound variants
  variant_ptr = variant_ptr->Union(compound_ptr);

  ExecEnv::log().info("Generated: {} Total variants for Genome: {}", variant_ptr->variantCount(), genome_name);

  return variant_ptr;

}


size_t kgl::VariantFactory::addGenomeSingleThreadVariant(std::shared_ptr<GenomeVariant> genome_variants,
                                             std::shared_ptr<const Variant> variant_ptr) const {

  genome_variants->addVariant(variant_ptr);
  return 1;

}



