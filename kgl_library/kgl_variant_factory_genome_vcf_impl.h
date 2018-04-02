//
// Created by kellerberrin on 27/02/18.
//

#ifndef KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H

#include "kgl_variant_factory_vcf_impl.h"
#include "kgl_variant_factory_vcf_cigar_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Multi threaded implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParseGenomeVCFImpl : public ParseCigarImpl {

public:

  ParseGenomeVCFImpl(const std::string &genome_name,
                     std::shared_ptr<PopulationVariant> pop_variant_ptr,
                     std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                     const std::string &vcf_file_name,
                     Phred_t variant_quality) : ParseCigarImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality),
                                                genome_name_(genome_name)  {}

  ~ParseGenomeVCFImpl() override = default;

  void readParseVCFImpl() override;

  void ProcessVCFRecord(const seqan::VcfRecord& vcf_record) override;

protected:

  const std::string genome_name_;

private:

  virtual bool parseVcfRecord(const std::string& genome_name,
                              const seqan::VcfRecord& record,
                              std::shared_ptr<const ContigFeatures> contig_ptr,
                              Phred_t variant_quality,
                              bool& quality_ok,
                              size_t& record_variants) = 0;


};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H
