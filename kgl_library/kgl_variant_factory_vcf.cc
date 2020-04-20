//
// Created by kellerberrin on 26/12/17.
//

#include "kel_exec_env.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_pf3k_impl.h"
#include "kgl_variant_factory_grch_impl.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto the implementation objects.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::VcfFactory::parseVCFVariants( std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                        std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                        const std::string &vcf_file_name,
                                        VCFParserEnum parser_type) {

  switch(parser_type) {

    case VCFParserEnum::GatkMultiGenome:
      return gatkMultiGenomeVCFVariants(vcf_population_ptr, genome_db_ptr, vcf_file_name);

    case VCFParserEnum::GRChNoGenome:
      return GRChNoGenomeVCFVariants(vcf_population_ptr, genome_db_ptr, vcf_file_name);

    case VCFParserEnum::NotImplemented:
      ExecEnv::log().error("Variant Parser not specified or not implemented, VCF file: {}", vcf_file_name);
      break;

  }

  return false;

}



bool kgl::VcfFactory::gatkMultiGenomeVCFVariants(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                                 std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                 const std::string &vcf_file_name) {

  Pf3kVCFImpl reader(vcf_population_ptr, genome_db_ptr, vcf_file_name);

  reader.readParseVCFImpl();

  return true;

}

bool kgl::VcfFactory::GRChNoGenomeVCFVariants(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                              const std::string &vcf_file_name) {

  GrchVCFImpl reader(vcf_population_ptr, genome_db_ptr, vcf_file_name);

  reader.readParseVCFImpl();

  return true;

}

