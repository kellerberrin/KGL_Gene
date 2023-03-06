//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kel_utility.h"
#include "kgl_variant_db_population.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_parse_info.h"



namespace kellerberrin::genome {   //  organization level namespace


class PfVCFImpl : public VCFReaderMT {

public:

  PfVCFImpl(const std::shared_ptr<PopulationDB>& vcf_population_ptr,
            const std::shared_ptr<const GenomeReference>& genome_db_ptr,
            const ContigAliasMap&,    // Chromosome aliasing is not used on Gatk (P.Falciparum) VCF files.
            const EvidenceInfoSet& evidence_map) : evidence_factory_(evidence_map),
                                                   unphased_population_ptr_(vcf_population_ptr),
                                                   genome_db_ptr_(genome_db_ptr) {}
  ~PfVCFImpl() override = default;

  void ProcessVCFRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr) override;

  void processVCFHeader(const VCFHeaderInfo& header_info) override;

  void readParseVCFImpl(const std::string &vcf_file_name);

private:

  EvidenceFactory evidence_factory_;

  // Processes the record in a try/catch block.
  void ParseRecord(std::unique_ptr<const VCFRecord> vcf_record_ptr);

  constexpr static const char GT_FIELD_SEPARATOR_CHAR_{'/'};
  constexpr static const char AD_FIELD_SEPARATOR_CHAR_{','};
  constexpr static const char* DIGITS_{"0123456789"};
  constexpr static const char* FLOAT_DIGITS_{"eE.+-0123456789"};
  constexpr static const char* UPSTREAM_ALLELE_{"*"};
  constexpr static const char FORMAT_SEPARATOR_ = ':';
  constexpr static const char* MISSING_VALUE_{"."};
  constexpr static const size_t DIPLOID_{2};


  // Progress counters.
  size_t variant_count_{0};
  constexpr static const size_t VARIANT_REPORT_INTERVAL_{10000};

  // This object is write accessed by multiple threads, it MUST BE mutex guarded for any access.
  const std::shared_ptr<PopulationDB> unphased_population_ptr_;   // Un-phased variants.
  const std::shared_ptr<const GenomeReference> genome_db_ptr_; // read access only.

  void setupPopulationStructure(const std::shared_ptr<const GenomeReference>& genome_db_ptr);

  [[nodiscard]] bool createAddVariant(const std::string& genome_name,
                                      const std::shared_ptr<const ContigReference>& contig_ptr,
                                      ContigOffset_t contig_offset,
                                      const std::string& identifier,
                                      const std::string& reference,
                                      const std::string& alternate,
                                      const VariantEvidence& evidence);

  bool addThreadSafeVariant(const std::shared_ptr<const Variant>& variant_ptr, const GenomeId_t& genome) const;

};



}   // end namespace



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
