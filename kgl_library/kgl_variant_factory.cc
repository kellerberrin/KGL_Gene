//
// Created by kellerberrin on 27/11/17.
//

#include "kel_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;



void kgl::VariantFactory::readVCFVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                          std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                          const std::string& variant_file_name) const {


  ExecEnv::log().info("Processing VCF file: {}", variant_file_name);

  if (not VcfFactory().readParseVCFVariants(vcf_population_ptr, genome_db_ptr, variant_file_name)) {

    ExecEnv::log().error("Problem parsing VCF file: {}", variant_file_name);

  }

}


