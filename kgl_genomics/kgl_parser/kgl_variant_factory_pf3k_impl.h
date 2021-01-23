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


class Pf3kVCFImpl : public VCFReaderMT {

public:

  Pf3kVCFImpl(const std::shared_ptr<PopulationDB> vcf_population_ptr,
              const std::shared_ptr<const GenomeReference> genome_db_ptr,
              const ContigAliasMap&,    // Chromosome aliasing is not used on Gatk (P.Falciparum) VCF files.
              const EvidenceInfoSet& evidence_map) : evidence_factory_(evidence_map),
                                                     unphased_population_ptr_(vcf_population_ptr),
                                                     genome_db_ptr_(genome_db_ptr) {}
  ~Pf3kVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo& header_info) override;

  void readParseVCFImpl(const std::string &vcf_file_name);

private:

  EvidenceFactory evidence_factory_;

  // Processes the record in a try/catch block.
  void ParseRecord(size_t vcf_record_count, const VcfRecord& record);

  constexpr static const char PL_CHECK_ZERO_ = '0';  // Check if the first PL character is zero, discard if true.
  constexpr static const char PL_CHECK_DOT_ = '.';  // Check if the first PL character is '.', discard if true.
  constexpr static const char* VQSLOD_INFO_FIELD_ = "VQSLOD";  // In the Info record.
  constexpr static const char* GT_FIELD_SEPARATOR_ = "/";
  constexpr static const char GT_FIELD_SEPARATOR_CHAR_ = '/';
  constexpr static const char* PL_FIELD_SEPARATOR_ = ",";
  constexpr static const char* AD_FIELD_SEPARATOR_ = ",";
  constexpr static const char AD_FIELD_SEPARATOR_CHAR_ = ',';
  constexpr static const char* DIGITS_ = "0123456789";
  constexpr static const char* FLOAT_DIGITS_ = "+-.e0123456789";
  constexpr static const char* UPSTREAM_ALLELE_ = "*";
  constexpr static const char* PASSED_FILTERS_{"PASS"};


// Quality constants.

  constexpr static const double MIN_DEPTH_ = 5;
  constexpr static const double MIN_QUALITY_ = 10.0;
  constexpr static const double MIN_VQSLOD_QUALITY_ = 5.0;
  constexpr static const double MIN_GQ_QUALITY_ = 10.0;

// Progress counters.

  size_t variant_count_{0};
  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;

  // This object is write accessed by multiple threads, it MUST BE mutex guarded for any access.
  const std::shared_ptr<PopulationDB> unphased_population_ptr_;   // Un-phased variants.
  const std::shared_ptr<const GenomeReference> genome_db_ptr_; // read access only.

  void setupPopulationStructure(const std::shared_ptr<const GenomeReference> genome_db_ptr);

  [[nodiscard]] bool createAddVariant(const std::string& genome_name,
                                      const std::shared_ptr<const ContigReference> contig_ptr,
                                      ContigOffset_t contig_offset,
                                      bool passed_filter,
                                      const std::string& identifier,
                                      const std::string& reference,
                                      const std::string& alternate,
                                      const VariantEvidence& evidence);

  bool addThreadSafeVariant(std::unique_ptr<const Variant>&&, GenomeId_t genome) const;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses the seqan::VcfRecord for each genotype.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


}   // end namespace



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
