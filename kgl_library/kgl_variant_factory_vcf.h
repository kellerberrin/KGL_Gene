//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_H
#define KGL_VARIANT_FACTORY_VCF_H

#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"
#include "kgl_logging.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class VcfFactory {

public:

  explicit VcfFactory();
  ~VcfFactory();

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeVariant> readParseFreeBayesVcf(const std::string &genome_name,
                                                       std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                       const std::string &vcf_file_name,
                                                       Phred_t variant_quality) const;

  // Functionality passed to the implmentation.
  std::shared_ptr<GenomeVariant> readParseGATKVcf(const std::string &genome_name,
                                                  std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                                  const std::string &vcf_file_name,
                                                  Phred_t variant_quality) const;


private:

  class ParseVCFImpl;           // Forward declaration of the VCF File reader implementation Base class
  class FreeBayesVCFImpl;       // Forward declaration of the Free Bayes VCF File reader implementation class
  std::unique_ptr<FreeBayesVCFImpl> fb_vcf_impl_ptr_;    // Read Free Bayes VCF file PIMPL
  class GATKVCFImpl;       // Forward declaration of the GATK VCF File reader implementation class
  std::unique_ptr<GATKVCFImpl> gatk_vcf_impl_ptr_;    // Read GATK VCF file PIMPL


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_VCF_H
