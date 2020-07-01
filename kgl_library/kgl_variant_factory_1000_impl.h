//
// Created by kellerberrin on 27/6/20.
//

#ifndef KGL_VARIANT_FACTORY_1000_IMPL_H
#define KGL_VARIANT_FACTORY_1000_IMPL_H



#include "kel_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_info.h"
#include "kgl_variant_factory_population.h"


namespace kellerberrin::genome {   //  organization level namespace


class Genome1000VCFImpl : public VCFReaderMT {

public:

  Genome1000VCFImpl(const std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                    const std::shared_ptr<const GenomeReference> genome_db_ptr,
                    const std::string &vcf_file_name,
                    const ContigAliasMap& contig_alias_map,
                    const EvidenceInfoSet& evidence_map) : VCFReaderMT(vcf_file_name),
                                                           evidence_factory_(evidence_map),
                                                           contig_alias_map_(contig_alias_map),
                                                           unphased_population_ptr_(vcf_population_ptr),
                                                           genome_db_ptr_(genome_db_ptr) {}
  ~Genome1000VCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo& header_info) override;

  void readParseVCFImpl();

private:

  EvidenceFactory evidence_factory_;
  ContigAliasMap contig_alias_map_;

  // Processes the record in a try/catch block.
  void ParseRecord(size_t vcf_record_count, const VcfRecord& record);

  // Progress counters.
  mutable size_t abstract_variant_count_{0};
  size_t actual_variant_count_{0};
  size_t variant_count_{0};

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{10000};

  constexpr static const size_t REFERENCE_VARIANT_INDEX_{0};
  constexpr static const char* REFERENCE_VARIANT_INDICATOR_{"."};
  constexpr static const char PHASE_MARKER_ {'|'};
  constexpr static const char MULIPLE_ALT_SEPARATOR_{','};
  constexpr static const char ABSTRACT_ALT_BRACKET_{'<'};

  const std::shared_ptr<UnphasedPopulation> unphased_population_ptr_;   // Un-phased variants.
  const std::shared_ptr<const GenomeReference> genome_db_ptr_; // read access only.

  // mutex to lock the UnphasedPopulation structure.
  mutable std::mutex add_variant_mutex_;

  bool addThreadSafeVariant(std::unique_ptr<const Variant>&&, const std::vector<GenomeId_t>& genome_vector) const;
// Calculates alternate indexes for the two phases (.first = A, .second = B).
  std::pair<size_t, size_t> alternateIndex(const std::string& genotype, const std::vector<std::string>& alt_vector) const;
// Adds variants to a vector of genomes.
  void addVariants( const std::map<size_t, std::vector<GenomeId_t>>& phase_map,
                    const ContigId_t& contig,
                    PhaseId_t phase,
                    ContigOffset_t offset,
                    const InfoDataEvidence& info_evidence_opt,
                    const std::string& reference,
                    const std::vector<std::string>& alt_vector,
                    size_t vcf_record_count);

};



}   // end namespace



#endif //KGL_VARIANT_FACTORY_1000_IMPL_H
