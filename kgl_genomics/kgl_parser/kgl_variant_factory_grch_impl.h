//
// Created by kellerberrin on 20/4/20.
//

#ifndef KGL_VARIANT_FACTORY_GRCH_IMPL_H
#define KGL_VARIANT_FACTORY_GRCH_IMPL_H




#include "kel_utility.h"
#include "kgl_variant_db_genome.h"
#include "kgl_data_file_impl.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_db_population.h"


namespace kellerberrin::genome {   //  organization level namespace




class GrchVCFImpl : public VCFReaderMT {

public:

  GrchVCFImpl(const std::shared_ptr<PopulationDB>& population_ptr,
              const std::shared_ptr<const GenomeReference>& genome_db_ptr,
              const ContigAliasMap& contig_alias_map,
              const EvidenceInfoSet& evidence_map) : unphased_population_ptr_(population_ptr),
                                                     genome_db_ptr_(genome_db_ptr),
                                                     contig_alias_map_(contig_alias_map),
                                                     evidence_factory_(evidence_map) {}

  ~GrchVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord &vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo &header_info) override;

  void readParseVCFImpl(const std::string &vcf_file_name);


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{100000};
  constexpr static const char MULIPLE_ALT_SEPARATOR_{','};
  constexpr static const char* PASSED_FILTERS_{"PASS"};

  const std::shared_ptr<PopulationDB> unphased_population_ptr_;   // Un-phased variants.
  std::shared_ptr<const GenomeReference> genome_db_ptr_;
  ContigAliasMap contig_alias_map_;
  EvidenceFactory evidence_factory_;

// Progress counters.
  size_t variant_count_{0};

  bool addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr, const GenomeId_t& genome) const;
// Filters the variants to be added to the population structure.
  virtual bool filterVariant(const std::shared_ptr<const Variant>&,  const GenomeId_t&) const { return true; }

};




class SNPdbVCFImpl : public GrchVCFImpl {

public:

  SNPdbVCFImpl( const std::shared_ptr<PopulationDB>& population_ptr,
                const std::shared_ptr<const GenomeReference>& genome_db_ptr,
                const ContigAliasMap& contig_alias_map,
                const EvidenceInfoSet& evidence_map) : GrchVCFImpl(population_ptr, genome_db_ptr, contig_alias_map, evidence_map) {}

  ~SNPdbVCFImpl() override = default;

private:

  // Filters the variants to be added to the population structure.
  bool filterVariant(const std::shared_ptr<const Variant>& variant_ptr,  const GenomeId_t& genome) const override;

};


} // namespace

#endif //KGL_KGL_VARIANT_FACTORY_GRCH_IMPL_H
