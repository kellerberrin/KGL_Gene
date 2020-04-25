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




std::shared_ptr<kgl::UnphasedPopulation>
kgl::VcfFactory::gatkMultiGenomeVCFVariants( const std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                             const std::string &vcf_file_name) {

  std::shared_ptr<UnphasedPopulation> vcf_population_ptr(std::make_shared<UnphasedPopulation>(genome_db_ptr->genomeId()));
  Pf3kVCFImpl reader(vcf_population_ptr, genome_db_ptr, vcf_file_name);

  reader.readParseVCFImpl();

  return vcf_population_ptr;

}

std::shared_ptr<kgl::UnphasedGenome>
kgl::VcfFactory::GRChNoGenomeVCFVariants( const std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                          const std::string &vcf_file_name,
                                          const ContigAliasMap& contig_alias_map) {

  std::shared_ptr<UnphasedGenome> vcf_genome_ptr(std::make_shared<UnphasedGenome>(genome_db_ptr->genomeId()));
  GrchVCFImpl reader(vcf_genome_ptr, genome_db_ptr, vcf_file_name, contig_alias_map);

  reader.readParseVCFImpl();

  return vcf_genome_ptr;

}

