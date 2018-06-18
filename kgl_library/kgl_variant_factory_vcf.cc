//
// Created by kellerberrin on 26/12/17.
//

#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_pf3k_impl.h"


namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VcfFactory() is a public facade class that passes the functionality onto the implementation objects.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




bool kgl::VcfFactory::readParseVCFVariants(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                                            std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                            const std::string &vcf_file_name,
                                            Phred_t variant_quality) const {

  Pf3kVCFImpl reader(vcf_population_ptr, genome_db_ptr, vcf_file_name, variant_quality);

  reader.readParseVCFImpl();

  return true;

}


