//
// Created by kellerberrin on 20/12/20.
//

#ifndef KGL_VARIANT_FACTORY_GNOMAD_IMPL_H
#define KGL_VARIANT_FACTORY_GNOMAD_IMPL_H


#include "kel_utility.h"
#include "kgl_variant_db_population.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_info.h"


namespace kellerberrin::genome {   //  organization level namespace


class GenomeGnomadVCFImpl : public VCFReaderMT {

public:

  GenomeGnomadVCFImpl(const std::shared_ptr<PopulationDB> vcf_population_ptr,
                      const std::shared_ptr<const GenomeReference> genome_db_ptr,
                      const ContigAliasMap &contig_alias_map,
                      const EvidenceInfoSet &evidence_map) : evidence_factory_(evidence_map),
                                                             contig_alias_map_(contig_alias_map),
                                                             population_ptr_(vcf_population_ptr),
                                                             genome_db_ptr_(genome_db_ptr) {}

  ~GenomeGnomadVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord &vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo &header_info) override;

  void readParseVCFImpl(const std::string &vcf_file_name);

private:

  EvidenceFactory evidence_factory_;
  ContigAliasMap contig_alias_map_;

  // Processes the record in a try/catch block.
  void ParseRecord(size_t vcf_record_count, const VcfRecord &record);

  // Progress counters.
  size_t actual_variant_count_{0};
  size_t variant_count_{0};

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{10000};

  constexpr static const size_t REFERENCE_VARIANT_INDEX_{0};
  constexpr static const char *REFERENCE_VARIANT_INDICATOR_{"."};
  constexpr static const char GT_SEPARATOR_{':'};
  constexpr static const char PHASE_MARKER_{'/'};
  constexpr static const size_t MINIMUM_GENOTYPE_SIZE_{3};
  constexpr static const char MULTIPLE_ALT_SEPARATOR_{','};
  constexpr static const char ABSTRACT_ALT_BRACKET_{'<'};
  constexpr static const char *PASSED_FILTERS_{"PASS"};

  const std::shared_ptr<PopulationDB> population_ptr_;   // Diploid phased variants.
  const std::shared_ptr<const GenomeReference> genome_db_ptr_; // read access only.

  // mutex to lock the structure for multiple thread access by parsers.
  mutable std::mutex add_variant_mutex_;

  bool addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr, const std::vector<GenomeId_t>& genome_vector) const;

  // Calculates alternate indexes for the two phases (.first = A, .second = B).
  std::pair<size_t, size_t> alternateIndex(const std::string &genotype, const std::vector<std::string> &alt_vector) const;

  // Adds variants to a vector of genomes.
  void addVariants(const std::map<size_t, std::vector<GenomeId_t>> &phase_map,
                   const ContigId_t &contig,
                   VariantPhase phase,
                   ContigOffset_t offset,
                   bool passedFilters,
                   const std::shared_ptr<const DataMemoryBlock>& info_evidence_ptr,
                   const std::string &reference,
                   const std::string &identifier,
                   const std::vector<std::string> &alt_vector,
                   size_t vcf_record_count);

};


} // namespace

#endif //KGL_VARIANT_FACTORY_GNOMAD_IMPL_H
