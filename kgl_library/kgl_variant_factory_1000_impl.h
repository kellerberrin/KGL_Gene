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
                                                           index_variants_(vcf_population_ptr),
                                                           population_ptr_(vcf_population_ptr),
                                                           genome_db_ptr_(genome_db_ptr) {}
  ~Genome1000VCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo& header_info) override;

  void readParseVCFImpl();

private:

  EvidenceFactory evidence_factory_;
  ContigAliasMap contig_alias_map_;
  IndexVariants index_variants_;   // Single threaded variant indexer.

  // Processes the record in a try/catch block.
  void ParseRecord(size_t vcf_record_count, const VcfRecord& record);

  // Progress counters.
  mutable size_t abstract_variant_count_{0};
  size_t variant_count_{0};

  constexpr static const size_t VARIANT_REPORT_INTERVAL_{10000};

  constexpr static const size_t REFERENCE_VARIANT_INDEX_{0};
  constexpr static const char* REFERENCE_VARIANT_INDICATOR_{"."};
  constexpr static const char PHASE_MARKER_ {'|'};
  constexpr static const char MULIPLE_ALT_SEPARATOR_{','};
  constexpr static const char ABSTRACT_ALT_BRACKET_{'<'};

  // This object is write accessed by multiple threads, it MUST BE mutex guarded for any access.
  const std::shared_ptr<UnphasedPopulation> population_ptr_;   // Phased variants.
  const std::shared_ptr<const GenomeReference> genome_db_ptr_; // read access only.

// Calculates alternate indexes for the two phases (.first = A, .second = B).
  std::pair<size_t, size_t> alternateIndex(const std::string& genotype, const std::vector<std::string>& alt_vector) const;

};



}   // end namespace



#endif //KGL_VARIANT_FACTORY_1000_IMPL_H
