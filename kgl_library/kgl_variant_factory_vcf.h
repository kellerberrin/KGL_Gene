//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_H
#define KGL_VARIANT_FACTORY_VCF_H

#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"
#include "kgl_variant_factory_vcf_phasing.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class VcfFactory  {

public:

  VcfFactory() = default;
  ~VcfFactory() = default;

  // Functionality passed to the implementation.
  bool readParseFreeBayesVcf(const std::string &genome_name,
                             std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                             const std::string &vcf_file_name,
                             Phred_t variant_quality) const;

  bool readParseGATKVcf(const std::string &genome_name,
                        std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                        std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                        const std::string &vcf_file_name,
                        Phred_t variant_quality) const;

  bool readParsePf3kVariants(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                             std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                             const std::string &vcf_file_name,
                             Phred_t variant_quality) const;

private:


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_VCF_H
