//
// Created by kellerberrin on 27/11/17.
//


#include "kgl_variant_single.h"
#include "kgl_variant_factory.h"
#include "kgl_variant_factory_compound.h"


namespace kgl = kellerberrin::genome;


// Generate all SNP variants.
std::shared_ptr<const kgl::GenomeVariant>
kgl::VariantFactory::createVariants(const std::string &genome_name,
                                    const std::shared_ptr<const ContigCountData> &count_data,
                                    const std::shared_ptr<const GenomeDatabase> &genome_db_ptr,
                                    NucleotideReadCount_t minimum_read_count,
                                    double minimum_proportion) const {


  ExecEnv::log().info("Generating SNP variants for Genome: {}", genome_name);
  // generate snp raw variants.
  std::shared_ptr<const GenomeVariant> snp_ptr = SNPFactory().createSNPs(genome_name,
                                                                         count_data,
                                                                         genome_db_ptr,
                                                                         minimum_read_count,
                                                                         minimum_proportion);

  ExecEnv::log().info("Generated: {} SNP variants for Genome: {}", snp_ptr->size(), genome_name);
  ExecEnv::log().info("Generating Compound variants for Genome: {}", genome_name);

  // Generate contiguous deletion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_delete_ptr = kgl::CompoundDeleteFactory().create(snp_ptr,
                                                                                                 genome_db_ptr);
  // Generate contiguous insertion variants.
  std::shared_ptr<const kgl::GenomeVariant> cmp_insert_ptr = kgl::CompoundInsertFactory().create(snp_ptr,
                                                                                                 genome_db_ptr);
  // Generate compound single codon variants
  std::shared_ptr<const kgl::GenomeVariant> cmp_snp_ptr = kgl::CompoundSNPFactory().create(snp_ptr,
                                                                                           genome_db_ptr);

  // combine the compound variants
  std::shared_ptr<const kgl::GenomeVariant> compound_ptr = cmp_delete_ptr->Union(cmp_insert_ptr);
  compound_ptr = compound_ptr->Union(cmp_snp_ptr);

  ExecEnv::log().info("Generated: {} Compound variants for Genome: {}", compound_ptr->size(), genome_name);
  ExecEnv::log().info("Combining Compound and SNP variants for Genome: {}", genome_name);

  // disaggregate the compound variants.
  std::shared_ptr<const kgl::GenomeVariant> disagg_ptr = kgl::CompoundDeleteFactory().disaggregate(compound_ptr,
                                                                                                   genome_db_ptr);

  // remove the disaggregated snps from the set of raw snp variants.
  snp_ptr = snp_ptr->Difference(disagg_ptr);

  // add in the compound variants
  snp_ptr = snp_ptr->Union(compound_ptr);

  ExecEnv::log().info("Generated: {} Total variants for Genome: {}", snp_ptr->size(), genome_name);

  return snp_ptr;

}
