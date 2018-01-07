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

    variant_ptr = createVcfVariants(genome_db_ptr,
                                    genome_name,
                                    variant_file_name,
                                    variant_quality);

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
kgl::VariantFactory::createVcfVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                       const std::string& genome_name,
                                       const std::string& vcf_file_name,
                                       Phred_t variant_quality) const {


  ExecEnv::log().info("Generating VCF file: {}  variants for Genome: {}", vcf_file_name, genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> single_variant_ptr = VcfFactory().readParseVcf(genome_name,
                                                                                      genome_db_ptr,
                                                                                      vcf_file_name,
                                                                                      variant_quality);

  ExecEnv::log().info("Generated: {} Single variants for Genome: {}", single_variant_ptr->size(), genome_name);

  std::shared_ptr<const GenomeVariant> variant_ptr = aggregateVariants(genome_db_ptr, genome_name, single_variant_ptr);

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
  // Generate compound single codon variants
  std::shared_ptr<const kgl::GenomeVariant> cmp_snp_ptr = kgl::CompoundSNPFactory().create(single_variant_ptr,
                                                                                           genome_db_ptr);

  // Generate compound non-coding insert variants.
  std::shared_ptr<const kgl::GenomeVariant> non_coding_insert_ptr = kgl::CompoundNonCodingInsertFactory().create(single_variant_ptr,
                                                                                                                 genome_db_ptr);
  // Generate compound non-coding insert variants.
  std::shared_ptr<const kgl::GenomeVariant> non_coding_delete_ptr = kgl::CompoundNonCodingDeleteFactory().create(single_variant_ptr,
                                                                                                                 genome_db_ptr);

  // combine the compound variants
  std::shared_ptr<const kgl::GenomeVariant> compound_ptr = cmp_delete_ptr->Union(cmp_insert_ptr);
  compound_ptr = compound_ptr->Union(cmp_snp_ptr);
  compound_ptr = compound_ptr->Union(non_coding_insert_ptr);
  compound_ptr = compound_ptr->Union(non_coding_delete_ptr);

  ExecEnv::log().info("Generated: {} Compound variants for Genome: {}", compound_ptr->size(), genome_name);
  ExecEnv::log().info("Combining Compound and SNP variants for Genome: {}", genome_name);

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


// This function will insert multiple variants for each CDS sequence within each gene.
size_t kgl::VariantFactory::addSingleVariant(std::shared_ptr<GenomeVariant> genome_single_variants,
                                             std::shared_ptr<Variant> variant_ptr) {

  // Annotate the variant with genome information.
  size_t variant_count = 0;
  GeneVector gene_vector;
  ContigOffset_t variant_offset = variant_ptr->contigOffset();
  if (variant_ptr->contig()->findGenes(variant_offset, gene_vector)) {

    for (const auto& gene_ptr : gene_vector) {

      std::shared_ptr<const CodingSequenceArray> sequence_array = kgl::GeneFeature::getCodingSequences(gene_ptr);
      if (sequence_array->empty()) {

        std::shared_ptr<Variant> intron_variant_ptr = variant_ptr->clone();
        intron_variant_ptr->defineIntron(gene_ptr); // intron
        genome_single_variants->addVariant(variant_ptr);
        ++variant_count;

      } else {

        for (const auto& sequence : sequence_array->getMap()) {

          if (sequence.second->isWithinCoding(variant_offset)) {

            // create a variant copy and annotate with a coding sequence.
            std::shared_ptr<Variant> coding_single_ptr = variant_ptr->clone();
            coding_single_ptr->defineCoding(sequence.second); // coding
            genome_single_variants->addVariant(coding_single_ptr);
            ++variant_count;

          } else {  // an intron for this sequence

            // create a variant copy and annotate with a gene.
            std::shared_ptr<Variant> intron_single_ptr = variant_ptr->clone();
            intron_single_ptr->defineIntron(gene_ptr); // intron
            genome_single_variants->addVariant(intron_single_ptr);
            ++variant_count;

          } // if valid sequence for offset

        } // for all sequences within a gene

      } // if gene has a valid sequence.

    } // for all genes.

  } else {

    // create a variant copy and tag as non-coding.
    variant_ptr->defineNonCoding(); // non coding
    genome_single_variants->addVariant(variant_ptr);
    ++variant_count;

  }

  return variant_count;

}

