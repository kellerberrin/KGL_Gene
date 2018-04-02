//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_FBVCF_IMPL_H
#define KGL_VARIANT_FACTORY_FBVCF_IMPL_H




#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_vcf_parse_impl.h"
#include "kgl_variant_factory_genome_vcf_impl.h"


#include <seqan/vcf_io.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF (freebayes) parser. Low-level implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class FreeBayesVCFImpl : public ParseGenomeVCFImpl {

public:

  FreeBayesVCFImpl(const std::string& genome_name,
                   std::shared_ptr<PopulationVariant> pop_variant_ptr,
                   std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                   const std::string& vcf_file_name,
                   Phred_t variant_quality) : ParseGenomeVCFImpl(genome_name,
                                                                 pop_variant_ptr,
                                                                 genome_db_ptr,
                                                                 vcf_file_name,
                                                                 variant_quality) {}

  ~FreeBayesVCFImpl() override = default;

private:


  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      std::shared_ptr<const ContigFeatures> contig_ptr,
                      Phred_t variant_quality,
                      bool& quality_ok,
                      size_t& variant_count) override;


};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_FBVCF_IMPL_H
