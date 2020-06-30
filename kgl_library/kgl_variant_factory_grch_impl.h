//
// Created by kellerberrin on 20/4/20.
//

#ifndef KGL_VARIANT_FACTORY_GRCH_IMPL_H
#define KGL_VARIANT_FACTORY_GRCH_IMPL_H




#include "kel_utility.h"
#include "kgl_variant_db_unphased.h"
#include "kgl_variant_file_impl.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_population.h"



namespace kellerberrin::genome {   //  organization level namespace




class GrchVCFImpl : public VCFReaderMT {

public:

  GrchVCFImpl(std::shared_ptr<UnphasedPopulation> population_ptr,
              std::shared_ptr<const GenomeReference> genome_db_ptr,
              const std::string &vcf_file_name,
              const ContigAliasMap& contig_alias_map,
              const EvidenceInfoSet& evidence_map) : VCFReaderMT(vcf_file_name),
                                                        unphased_population_ptr_(std::move(population_ptr)),
                                                        genome_db_ptr_(std::move(genome_db_ptr)),
                                                        contig_alias_map_(contig_alias_map),
                                                        evidence_factory_(evidence_map) {}

  ~GrchVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord &vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo &header_info) override;

  void readParseVCFImpl();


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{100000};
  constexpr static const char MULIPLE_ALT_SEPARATOR_{','};

  const std::shared_ptr<UnphasedPopulation> unphased_population_ptr_;   // Un-phased variants.
  std::shared_ptr<const GenomeReference> genome_db_ptr_;
  ContigAliasMap contig_alias_map_;
  EvidenceFactory evidence_factory_;

// Progress counters.
  size_t variant_count_{0};
  // mutex to lock the UnphasedPopulation structure.
  mutable std::mutex add_variant_mutex_;

  bool addThreadSafeVariant(std::unique_ptr<const Variant>&&, GenomeId_t genome) const;

};

} // namespace

#endif //KGL_KGL_VARIANT_FACTORY_GRCH_IMPL_H
