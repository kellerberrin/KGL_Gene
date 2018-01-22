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


namespace kgl = kellerberrin::genome;


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createVariants(std::shared_ptr<const kgl::GenomeDatabase> genome_db_ptr,
                                    const std::string& genome_name,
                                    const std::string& variant_file_name,
                                    bool vcf_is_gatk,
                                    Phred_t read_quality,
                                    Phred_t variant_quality,
                                    NucleotideReadCount_t min_read_count,
                                    double min_proportion) const {



  std::string file_ext = kgl::Utility::fileExtension(variant_file_name);
  std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::toupper); // convert to UC for robust comparison

  std::shared_ptr<const kgl::GenomeVariant> variant_ptr;

  if (file_ext == SAM_FILE_EXTENSTION_) {

    variant_ptr = createSamVariants(genome_db_ptr,
                                    genome_name,
                                    variant_file_name,
                                    read_quality,
                                    variant_quality,
                                    min_read_count,
                                    min_proportion);

  } else if (file_ext == BAM_FILE_EXTENSTION_) {

    variant_ptr = createBamVariants(genome_db_ptr,
                                    genome_name,
                                    variant_file_name,
                                    read_quality,
                                    variant_quality,
                                    min_read_count,
                                    min_proportion);

  } else if (file_ext == VCF_FILE_EXTENSTION_) {

    if (vcf_is_gatk) {


      variant_ptr = createGATKVcfVariants(genome_db_ptr, genome_name, variant_file_name, variant_quality);


    } else { // Check if the prefix of the file is "gatk" (case insenstive).

      std::string file_name = kgl::Utility::fileName(variant_file_name);
      std::transform(file_name.begin(), file_name.end(), file_name.begin(), ::toupper); // convert to UC for robust comparison
      std::size_t gtk_prefix_length = std::strlen(GATK_FILE_PREFIX_);
      std::string file_name_prefix = file_name.substr(0, gtk_prefix_length);

      ExecEnv::log().info("VCF file name prefix: {}", file_name_prefix);

      if (file_name_prefix == GATK_FILE_PREFIX_) {

        variant_ptr = createGATKVcfVariants(genome_db_ptr, genome_name, variant_file_name, variant_quality);

      } else {

        variant_ptr = createFreeBayesVcfVariants(genome_db_ptr, genome_name, variant_file_name, variant_quality);

      }


    }



  } else {

    ExecEnv::log().critical("Unsupported file type: '{}' for variant calling. Must be SAM ('.sam'), BAM ('.bam') or freebayes VCF ('.vcf')", file_ext);

  }


  return variant_ptr;

}



std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createSamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::string& sam_file_name,
                                       Phred_t read_quality,
                                       Phred_t variant_quality,
                                       NucleotideReadCount_t min_read_count,
                                       double min_proportion) const {

  // Read in the SAM file.
  std::shared_ptr<const kgl::ContigCountData> count_data_ptr = kgl::SamCountReader().readSAMFile(genome_db_ptr,
                                                                                                 sam_file_name,
                                                                                                 read_quality);

  ExecEnv::log().info("Generating SAM file: {}  variants for Genome: {}", sam_file_name, genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> single_variant_ptr = SingleFactory().createSingleVariants(genome_name,
                                                                                                 count_data_ptr,
                                                                                                 genome_db_ptr,
                                                                                                 variant_quality,
                                                                                                 min_read_count,
                                                                                                 min_proportion);

  ExecEnv::log().info("Generated: {} Single variants for Genome: {}", single_variant_ptr->size(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, single_variant_ptr);

  return variant_ptr;

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createFreeBayesVcfVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                const std::string &genome_name,
                                                const std::string &vcf_file_name,
                                                Phred_t variant_quality) const {


  ExecEnv::log().info("Generating VCF file: {}  Freebayes variants for Genome: {}", vcf_file_name, genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> fb_variant_ptr = VcfFactory().readParseFreeBayesVcf(genome_name,
                                                                                           genome_db_ptr,
                                                                                           vcf_file_name,
                                                                                           variant_quality);

  ExecEnv::log().info("Generated: {} Freebayes variants for Genome: {}", fb_variant_ptr->size(), genome_name);

  // Do we still need this step?
  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, fb_variant_ptr);

  return variant_ptr;

}


std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createGATKVcfVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                           const std::string &genome_name,
                                           const std::string &vcf_file_name,
                                           Phred_t variant_quality) const {


  ExecEnv::log().info("Generating VCF file: {}  GATK variants for Genome: {}", vcf_file_name, genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> gatk_variant_ptr = VcfFactory().readParseGATKVcf(genome_name,
                                                                                        genome_db_ptr,
                                                                                        vcf_file_name,
                                                                                        variant_quality);

  ExecEnv::log().info("Generated: {} GATK variants for Genome: {}", gatk_variant_ptr->size(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, gatk_variant_ptr);

  return variant_ptr;

}



std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createBamVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::string& vcf_file_name,
                                       Phred_t read_quality,
                                       Phred_t variant_quality,
                                       NucleotideReadCount_t min_read_count,
                                       double min_proportion) const {


  ExecEnv::log().info("Generating BAM file: {}  variants for Genome: {}", vcf_file_name, genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> single_variant_ptr = BamFactory().readParseBam(genome_name,
                                                                                      genome_db_ptr,
                                                                                      vcf_file_name,
                                                                                      read_quality,
                                                                                      variant_quality,
                                                                                      min_read_count,
                                                                                      min_proportion);

  ExecEnv::log().info("Generated: {} Single variants for Genome: {}", single_variant_ptr->size(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, single_variant_ptr);

  return variant_ptr;

}


// Generate compound variants.
std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::aggregateVariants(const std::shared_ptr<const GenomeDatabase>& genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::shared_ptr<const kgl::GenomeVariant>& single_variant_ptr) const {


  ExecEnv::log().info("Generating Compound variants for Genome: {}", genome_name);

  // Generate contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_delete_ptr = kgl::CompoundDeleteFactory().create(single_variant_ptr,
                                                                                                 genome_db_ptr);
  // Generate contiguous insertion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_insert_ptr = kgl::CompoundInsertFactory().create(single_variant_ptr,
                                                                                                 genome_db_ptr);

  // combine the compound variants
  std::shared_ptr<const kgl::GenomeVariant> compound_ptr = cmp_delete_ptr->Union(cmp_insert_ptr);

  ExecEnv::log().info("Generated: {} Compound variants for Genome: {}", compound_ptr->size(), genome_name);
  ExecEnv::log().info("Combining Compound and Single variants for Genome: {}", genome_name);

  // disaggregate the compound variants.
  std::shared_ptr<const kgl::GenomeVariant> disagg_ptr = kgl::CompoundDeleteFactory().disaggregate(compound_ptr,
                                                                                                   genome_db_ptr);

  // remove the disaggregated snps from the set of raw snp variants.
  std::shared_ptr<const kgl::GenomeVariant> variant_ptr = single_variant_ptr->Difference(disagg_ptr);

  // add in the compound variants
  variant_ptr = variant_ptr->Union(compound_ptr);

  ExecEnv::log().info("Generated: {} Total variants for Genome: {}", variant_ptr->size(), genome_name);

  return variant_ptr;

}


size_t kgl::VariantFactory::addGenomeVariant(std::shared_ptr<GenomeVariant> genome_variants,
                                             std::shared_ptr<const Variant> variant_ptr) {

  genome_variants->addVariant(variant_ptr);
  return 1;

}

