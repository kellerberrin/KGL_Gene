//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Pf3kVCFImpl : public ParseVCFImpl {

public:

  Pf3kVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name,
              Phred_t variant_quality) : ParseVCFImpl(pop_variant_ptr, genome_db_ptr, vcf_file_name, variant_quality) {}
  ~Pf3kVCFImpl() = default;

  void ProcessVCFRecord(const seqan::VcfRecord& record_ptr) override;

private:

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 1000;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
