//
// Created by kellerberrin on 22/11/17.
//

#ifndef KGL_VARIANT_FACTORY_H
#define KGL_VARIANT_FACTORY_H



#include "kgl_variant_db.h"
#include "kgl_genome_db.h"
#include "kgl_variant.h"
#include "kgl_variant_factory_vcf_impl.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The top level variant factory object
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class VariantFactory {

public:

  explicit VariantFactory() = default;
  virtual ~VariantFactory() = default;

  void readVCFVariants(std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                       std::shared_ptr<UnphasedPopulation> vcf_population_ptr,
                       const std::string& variant_file_name) const;


private:

  constexpr static const char* VCF_FILE_EXTENSTION_ = ".VCF";

};



}   // namespace genome
}   // namespace kellerberrin



#endif //KGL_VARIANT_FACTORY_H
