//
// Created by kellerberrin on 27/11/17.
//

#include "kel_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory.h"
#include "kgl_filter.h"


namespace kgl = kellerberrin::genome;



std::shared_ptr<kgl::UnphasedPopulation> kgl::VariantFactory::readVCFVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                                              const std::string& variant_file_name,
                                                                              VCFParserEnum parser_type) {


  ExecEnv::log().info("Processing VCF file: {}", variant_file_name);

  std::shared_ptr<UnphasedPopulation> vcf_population_ptr(std::make_shared<UnphasedPopulation>());
  if (not VcfFactory::parseVCFVariants(vcf_population_ptr, genome_db_ptr, variant_file_name, parser_type)) {

    ExecEnv::log().error("Problem parsing VCF file: {}", variant_file_name);

  }

  return vcf_population_ptr;

}


