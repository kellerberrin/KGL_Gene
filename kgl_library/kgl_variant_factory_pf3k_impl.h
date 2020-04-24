//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kel_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_readvcf_impl.h"
#include "kgl_variant_factory_record_vcf_impl.h"



namespace kellerberrin::genome {   //  organization level namespace


class Pf3kVCFImpl : public VCFReaderMT {

public:

  Pf3kVCFImpl(const std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
              const std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name) : VCFReaderMT(vcf_file_name),
                                                  unphased_population_ptr_(vcf_population_ptr),
                                                  genome_db_ptr_(genome_db_ptr) {}
  ~Pf3kVCFImpl() override = default;

  void ProcessVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record) override;

  void processVCFHeader(const VcfHeaderInfo& header_info) override;

  void readParseVCFImpl();

private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;

  // Processes the record in a try/catch block.
  void TryVCFRecord(size_t vcf_record_count, const VcfRecord& vcf_record);
  void ParseRecord(size_t vcf_record_count, const VcfRecord& record);

  constexpr static const char PL_CHECK_ZERO_ = '0';  // Check if the first PL character is zero, discard if true.
  constexpr static const char PL_CHECK_DOT_ = '.';  // Check if the first PL character is '.', discard if true.
  constexpr static const char* VQSLOD_INFO_FIELD_ = "VQSLOD";  // In the Info record.
  constexpr static const char* GT_FIELD_SEPARATOR_ = "/";
  constexpr static const char* PL_FIELD_SEPARATOR_ = ",";
  constexpr static const char* AD_FIELD_SEPARATOR_ = ",";
  constexpr static const char* DIGITS_ = "0123456789";
  constexpr static const char* FLOAT_DIGITS_ = "+-.e0123456789";
  constexpr static const char* UPSTREAM_ALLELE_ = "*";


// Quality constants.

  constexpr static const double MIN_DEPTH_ = 5;
  constexpr static const double MIN_QUALITY_ = 10.0;
  constexpr static const double MIN_VQSLOD_QUALITY_ = 5.0;
  constexpr static const double MIN_GQ_QUALITY_ = 10.0;

// Progress counters.

  size_t vcf_variant_count_{0};
  std::atomic<uint64_t> record_count_{0};
  size_t variant_count_{0};

// Create and add variants to population pointer (thread safe).

  mutable std::mutex add_variant_mutex_;  // mutex to lock the UnphasedPopulation structure when inserting variants.
  // This object is write accessed by multiple threads, it MUST BE mutex guarded for any access.
  const std::shared_ptr<UnphasedPopulation> unphased_population_ptr_;   // Un-phased variants.
  const std::shared_ptr<const GenomeDatabase> genome_db_ptr_; // read access only.

  [[nodiscard]] bool addThreadSafeGenomeVariant(const std::shared_ptr<const Variant>& variant_ptr);
  void setupPopulationStructure(const std::shared_ptr<const GenomeDatabase> genome_db_ptr);
  [[nodiscard]] bool addThreadSafeVariant(std::shared_ptr<const Variant>& variant_ptr);

  [[nodiscard]] bool createAddVariant(const std::string& genome_name,
                                      const std::shared_ptr<const ContigFeatures> contig_ptr,
                                      ContigOffset_t contig_offset,
                                      const std::string& reference,
                                      const std::string& alternate,
                                      const std::shared_ptr<const VariantEvidence> evidence_ptr);

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses the seqan::VcfRecord for each genotype.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


}   // end namespace



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
