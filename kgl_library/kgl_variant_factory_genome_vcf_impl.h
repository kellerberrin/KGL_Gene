//
// Created by kellerberrin on 27/02/18.
//

#ifndef KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H
#define KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H

#include "kgl_variant_factory_vcf_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// VCF parser. Multi threaded implementation.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



class ParseGenomeVCFImpl : public ParseVCFImpl {

public:

  ParseGenomeVCFImpl(const std::string &genome_name,
                     std::shared_ptr<PopulationVariant> pop_variant_ptr,
                     std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                     const std::string &vcf_file_name,
                     Phred_t variant_quality) : ParseVCFImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality),
                                                genome_name_(genome_name)  {

    genome_single_variants_ = GenomeVariant::emptyGenomeVariant(genome_name_, genome_db_ptr_);

  }

  ~ParseGenomeVCFImpl() override = default;

  void readParseVCFImpl() override;

  void ProcessVCFRecord(const seqan::VcfRecord& record_ptr) override;

protected:

  const std::string genome_name_;
  std::shared_ptr<GenomeVariant> genome_single_variants_;

private:

  virtual bool parseVcfRecord(const std::string& genome_name,
                              const seqan::VcfRecord& record,
                              std::shared_ptr<const ContigFeatures> contig_ptr,
                              std::shared_ptr<GenomeVariant> genome_variants,
                              Phred_t variant_quality,
                              bool& quality_ok,
                              size_t& record_variants) const = 0;

  void addGenome(std::shared_ptr<const GenomeVariant> genome_variant_ptr,
                 std::shared_ptr<PopulationVariant> pop_variant_ptr,
                 Phred_t read_quality) const;

};


}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_VARIANT_FACTORY_GENOME_VCF_IMPL_H
