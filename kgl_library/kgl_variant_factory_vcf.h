//
// Created by kellerberrin on 26/12/17.
//

#ifndef KGL_VARIANT_FACTORY_VCF_H
#define KGL_VARIANT_FACTORY_VCF_H

#include <memory>
#include <string>
#include <map>
#include "kgl_variant_db.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace



class VcfFactory  {

public:

  VcfFactory() = default;
  ~VcfFactory() = default;

  bool readParseVCFVariants(std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                            std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                            const std::string &vcf_file_name) const;

private:


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_VARIANT_FACTORY_VCF_H
