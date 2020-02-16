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

  std::string file_ext = Utility::fileExtension(variant_file_name);
  std::transform(file_ext.begin(), file_ext.end(), file_ext.begin(), ::toupper); // convert to UC for robust comparison

  if (file_ext == VCF_FILE_EXTENSTION_) {

    ExecEnv::log().info("Processing VCF file: {}", variant_file_name);
    VcfFactory().readParseVCFVariants(vcf_population_ptr, genome_db_ptr, variant_file_name);

  } else {

    ExecEnv::log().error("Invalid file name: {}", variant_file_name);
    ExecEnv::log().critical("Unsupported file type: '{}' for variant calling. Must be VCF ('.vcf')", file_ext);

  }

}


