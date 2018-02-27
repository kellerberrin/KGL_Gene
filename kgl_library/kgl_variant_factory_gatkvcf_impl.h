//
// Created by kellerberrin on 22/01/18.
//

#ifndef KGL_VARIANT_FACTORY_GATKVCF_IMPL_H
#define KGL_VARIANT_FACTORY_GATKVCF_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_genome_vcf_impl.h"

#include <seqan/vcf_io.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class GATKVCFImpl : public ParseGenomeVCFImpl {

public:

  GATKVCFImpl(const std::string& genome_name,
              std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string& vcf_file_name,
              Phred_t variant_quality) : ParseGenomeVCFImpl(genome_name,
                                                            pop_variant_ptr,
                                                            genome_db_ptr,
                                                            vcf_file_name,
                                                            variant_quality) {}

  ~GATKVCFImpl() override = default;

private:


  bool parseVcfRecord(const std::string& genome_name,
                      const seqan::VcfRecord& record,
                      std::shared_ptr<const ContigFeatures> contig_ptr,
                      std::shared_ptr<GenomeVariant> genome_variants,
                      Phred_t variant_quality,
                      bool& quality_ok,
                      size_t& record_variants) const override;

  bool parseSNP(const std::string& variant_source,
                std::shared_ptr<const ContigFeatures> contig_ptr,
                std::shared_ptr<GenomeVariant> genome_variants,
                Phred_t quality,
                const std::string& info,
                const std::string& reference,
                const std::string& alternate,
                ContigOffset_t contig_offset,
                size_t& variant_count) const;

  bool parseInsert(const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& variant_count) const;

  bool parseDelete(const std::string& variant_source,
                   std::shared_ptr<const ContigFeatures> contig_ptr,
                   std::shared_ptr<GenomeVariant> genome_variants,
                   Phred_t quality,
                   const std::string& info,
                   const std::string& reference,
                   const std::string& alternate,
                   ContigOffset_t contig_offset,
                   size_t& variant_count) const;



};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_GATKVCF_IMPL_H
