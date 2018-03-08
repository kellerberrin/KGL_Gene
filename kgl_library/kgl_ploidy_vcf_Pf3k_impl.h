//
// Created by kellerberrin on 7/03/18.
//

#ifndef KGL_PLOIDY_VCF_PF3K_IMPL_H
#define KGL_PLOIDY_VCF_PF3K_IMPL_H

#include "kgl_ploidy_analysis.h"
#include "kgl_genome_db.h"
#include "kgl_variant_db.h"
#include "kgl_variant_factory_record_vcf_impl.h"

#include <seqan/vcf_io.h>


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace

class ProcessPloidy {


public:

  ProcessPloidy(const std::vector<std::string>& genome_names,
                std::shared_ptr<const GenomeDatabase> genome_db_ptr) : genome_names_(genome_names),
                                                                       genome_db_ptr_(genome_db_ptr) {

    diploid_genotypes_.generateGenotypeVector(MAX_GENOTYPES_);
    ploidy_count_ = 0;

  }
  ~ProcessPloidy() = default;

  void PloidyVCFRecord(const seqan::VcfRecord& vcf_record,
                       std::shared_ptr<PloidyAnalysis> ploidy_ptr,
                       const ContigId_t& contig_id);

private:

  constexpr static const char PL_CHECK_ZERO_ = '0';  // Check if the first PL character is zero, discard if true.
  constexpr static const char PL_CHECK_DOT_ = '.';  // Check if the first PL character is '.', discard if true.
  constexpr static const char* VQSLOD_INFO_FIELD_ = "VQSLOD";  // In the Info record.
  constexpr static const char* GT_FIELD_SEPARATOR_ = "/";
  constexpr static const char* PL_FIELD_SEPARATOR_ = ",";
  constexpr static const char* AD_FIELD_SEPARATOR_ = ",";
  constexpr static const char* DIGITS_ = "0123456789";
  constexpr static const char* FLOAT_DIGITS_ = "+-.e0123456789";


// Quality constants.

  constexpr static const double MIN_DEPTH_ = 20;
  constexpr static const double MIN_QUALITY_ = 100.0;
  constexpr static const double MIN_VQSLOD_QUALITY_ = 0.0;
  constexpr static const double MIN_GQ_QUALITY_ = 0.0;
  constexpr static const double HQ_GQ_PLOIDY_ = 20;

// Quality counters.

  std::atomic<uint64_t> quality_failed_{0};
  std::atomic<uint64_t> vqslod_failed_{0};
  std::atomic<uint64_t> GQ_failed_{0};

  constexpr static const size_t MAX_GENOTYPES_ = 10; // maximum number of alleles per VCF record.=
  size_t ploidy_count_;
  DiploidGenotypes diploid_genotypes_;

  const std::vector<std::string>& genome_names_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;

  mutable std::mutex mutex_;

};



}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_PLOIDY_VCF_PF3K_IMPL_H
