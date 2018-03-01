//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_record_vcf_impl.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Pf3kVCFImpl : public ParseVCFImpl {

public:

  Pf3kVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name,
              Phred_t variant_quality) : ParseVCFImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality) {

    diploid_genotypes_.generateGenotypeVector(MAX_GENOTYPES_);

  }
  ~Pf3kVCFImpl() = default;

  void ProcessVCFRecord(const seqan::VcfRecord& vcf_record) override;

private:

  constexpr static const size_t MAX_GENOTYPES_ = 10; // maximum number of alleles per VCF record.
  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 1000;
  constexpr static const char PL_CHECK_ZERO_ = '0';  // Check if the first PL character is zero, discard if true.

  DiploidGenotypes diploid_genotypes_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses the seqan::VcfRecord for each genotype.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
