//
// Created by kellerberrin on 25/02/18.
//

#ifndef KGL_VARIANT_FACTORY_PF3K_IMPL_H
#define KGL_VARIANT_FACTORY_PF3K_IMPL_H



#include "kgl_utility.h"
#include "kgl_variant_factory_vcf.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class Pf3kVCFImpl  {

public:

  Pf3kVCFImpl(std::shared_ptr<PopulationVariant> pop_variant_ptr,
              std::shared_ptr<const GenomeDatabase> genome_db_ptr,
              const std::string &vcf_file_name,
              Phred_t variant_quality) : pop_variant_ptr_(pop_variant_ptr),
                                         genome_db_ptr_(genome_db_ptr),
                                         vcf_file_name_(vcf_file_name),
                                         variant_quality_(variant_quality) {

    vcf_variant_count_ = 0;

  }
  ~Pf3kVCFImpl() = default;

  bool readParsePf3kVariants();

  void ProcessVCFRecord(const seqan::VcfRecord& record_ptr);

private:


  std::shared_ptr<PopulationVariant> pop_variant_ptr_;
  std::shared_ptr<const GenomeDatabase> genome_db_ptr_;
  std::string vcf_file_name_;
  Phred_t variant_quality_;
  size_t vcf_variant_count_;

  constexpr static const size_t VARIANT_REPORT_INTERVAL_ = 1000;

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_PF3K_IMPL_H
