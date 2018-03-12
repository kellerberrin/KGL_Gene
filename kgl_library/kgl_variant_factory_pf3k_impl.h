//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_variant_factory_vcf_cigar_impl.h"
#include "kgl_ploidy_vcf_Pf3k_impl.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Pf3kVCFImpl : public ParseCigarImpl {

public:

  Pf3kVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name,
              Phred_t variant_quality) : ParseCigarImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality),
                                         process_ploidy_(getGenomeNames(), genome_db_ptr) {

  }
  ~Pf3kVCFImpl() = default;

  void ProcessVCFRecord(const seqan::VcfRecord& vcf_record) override;


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;
  ProcessPloidy process_ploidy_;

  // Processes the record in a try/catch block.
  void TryVCFRecord(const seqan::VcfRecord& vcf_record);
  void ParseRecord(const seqan::VcfRecord& vcf_record, const ContigId_t& contig_id);

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

  constexpr static const double MIN_DEPTH_ = 20;
  constexpr static const double MIN_QUALITY_ = 100.0;
  constexpr static const double MIN_VQSLOD_QUALITY_ = 0.0;
  constexpr static const double MIN_GQ_QUALITY_ = 20.0;

  std::atomic<uint64_t> record_count_{0};
  size_t variant_count_{0};


};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses the seqan::VcfRecord for each genotype.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
