//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"
#include "kgl_variant_factory_record_vcf_impl.h"
#include "kgl_ploidy_vcf_Pf3k_impl.h"



namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Pf3kVCFImpl : public ParseVCFImpl {

public:

  Pf3kVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name,
              Phred_t variant_quality) : ParseVCFImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality),
                                         process_ploidy_(getGenomeNames(), genome_db_ptr) {

  }
  ~Pf3kVCFImpl() = default;

  void ProcessVCFRecord(const seqan::VcfRecord& vcf_record) override;
  // Processes the record in a try/catch block.
  void TryVCFRecord(const seqan::VcfRecord& vcf_record);


private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 10000;
  ProcessPloidy process_ploidy_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Parses the seqan::VcfRecord for each genotype.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////


}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
